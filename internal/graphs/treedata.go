// Package containing all structs and functions related to graph-like data
// structures used in CAMUS such as quartets, networks, and trees
package graphs

import (
	"github.com/bits-and-blooms/bitset"
	"github.com/evolbioinfo/gotree/tree"
)

// Expanded tree struct containing necessary preprocessed data
type TreeData struct {
	tree.Tree
	Children       [][]*tree.Node      // Children for each node
	IdToNodes      []*tree.Node        // Mapping between id and node pointer
	quartetSet     [][]Quartet         // Quartets relevant for each subtree
	quartetCounts  *map[Quartet]uint32 // Count of each unique quartet topology
	Depths         []int               // Distance from all nodes to the root
	NumLeavesBelow []uint64            // Number of leaves below node
	NLeaves        int                 // Number of leaves
	leafsets       []*bitset.BitSet    // Leaves under each node
	lca            [][]int             // LCA for each pair of node id
	tipIndexMap    map[uint16]int      // Tip index to node id map
}

// Preprocess tree data and makes TreeData struct. Pass nil for qCounts if you
// don't need quartets.
func MakeTreeData(tre *tree.Tree, qCounts map[Quartet]uint32) *TreeData {
	children := children(tre)
	below := countLeavesBelow(tre, children)
	leafsets := calcLeafset(tre, children)
	lca := calcLCAs(tre, children)
	depths := calcDepths(tre)
	idMap := mapIdToNodes(tre)
	var qSets [][]Quartet
	if qCounts != nil {
		qSets = mapQuartetsToVertices(tre, qCounts, leafsets)
	}
	tipIndexMap := makeTipIndexMap(tre)
	return &TreeData{Tree: *tre,
		Children:       children,
		lca:            lca,
		leafsets:       leafsets,
		IdToNodes:      idMap,
		Depths:         depths,
		NumLeavesBelow: below,
		quartetSet:     qSets,
		quartetCounts:  &qCounts,
		tipIndexMap:    tipIndexMap,
		NLeaves:        len(tre.AllTipNames()),
	}
}

// Create mapping from id to node pointer
func mapIdToNodes(tre *tree.Tree) []*tree.Node {
	idMap := make([]*tree.Node, len(tre.Nodes()))
	tre.PostOrder(func(cur, prev *tree.Node, e *tree.Edge) (keep bool) {
		idMap[cur.Id()] = cur
		return true
	})
	return idMap
}

// Verify that tree still has the same root, and thus the data is still
// applicable
func (td *TreeData) Verify() {
	root := td.Root()
	if root != td.Root() {
		panic("TreeData root is wrong!")
	}
}

// Calculate children for each node for quick access (as gotree's Tree only
// stores neighbors)
func children(tre *tree.Tree) [][]*tree.Node {
	nNodes := len(tre.Nodes())
	children := make([][]*tree.Node, nNodes)
	tre.PostOrder(func(cur, prev *tree.Node, e *tree.Edge) (keep bool) {
		children[cur.Id()] = make([]*tree.Node, 2)
		if cur.Tip() {
			children[cur.Id()] = []*tree.Node{nil, nil}
		} else {
			children[cur.Id()] = GetChildren(cur)
		}
		return true
	})
	return children
}

// Get children of node
func GetChildren(node *tree.Node) []*tree.Node {
	children := make([]*tree.Node, 0)
	p, err := node.Parent()
	if err != nil && err.Error() == "The node has more than one parent" {
		panic(err)
	}
	i := 0
	for _, u := range node.Neigh() {
		if u != p {
			children = append(children, u)
			i++
		}
	}
	return children
}

// Calculates the leafset for every node
func calcLeafset(tre *tree.Tree, children [][]*tree.Node) []*bitset.BitSet {
	nLeaves, err := tre.NbTips()
	if err != nil {
		panic(err)
	}
	nNodes := len(tre.Nodes())
	leafset := make([]*bitset.BitSet, nNodes)
	tre.PostOrder(func(cur, prev *tree.Node, e *tree.Edge) (keep bool) {
		leafset[cur.Id()] = bitset.New(uint(nLeaves))
		if cur.Tip() {
			leafset[cur.Id()].Set(uint(cur.TipIndex()))
		} else {
			leafset[cur.Id()] = leafset[children[cur.Id()][0].Id()].Clone()
			for i := range len(children[cur.Id()]) - 1 {
				leafset[cur.Id()].InPlaceUnion(leafset[children[cur.Id()][i+1].Id()])
			}
		}
		return true
	})
	return leafset
}

// Calculates the LCA for every pair of nodes
func calcLCAs(tre *tree.Tree, children [][]*tree.Node) [][]int {
	nNodes := len(tre.Nodes())
	lca, below := make([][]int, nNodes), make([][]bool, nNodes) // below[i][j] = true means node j is below node i
	for i := range nNodes {
		lca[i] = make([]int, nNodes)
		below[i] = make([]bool, nNodes)
	}
	tre.PostOrder(func(cur, prev *tree.Node, e *tree.Edge) (keep bool) {
		below[cur.Id()][cur.Id()] = true
		lca[cur.Id()][cur.Id()] = cur.Id()
		if !cur.Tip() {
			for i := range nNodes {
				for _, child := range children[cur.Id()] {
					below[cur.Id()][i] = below[cur.Id()][i] || below[child.Id()][i]
				}
			}
			for c1 := range children[cur.Id()] {
				for c2 := c1 + 1; c2 < len(children[cur.Id()]); c2++ {
					for i := range nNodes {
						for j := range nNodes {
							childId1 := children[cur.Id()][c1].Id()
							childId2 := children[cur.Id()][c2].Id()
							if below[childId1][i] && below[childId2][j] ||
								i == cur.Id() && below[childId2][j] ||
								i == cur.Id() && below[childId1][j] {
								lca[i][j] = cur.Id()
								lca[j][i] = cur.Id()
							}
						}
					}
				}
			}
		}
		return true
	})
	return lca
}

// Calculate depths for all nodes in tree (slice index = node id)
func calcDepths(tre *tree.Tree) []int {
	depths := make([]int, len(tre.Nodes()))
	tre.PreOrder(func(cur, prev *tree.Node, e *tree.Edge) (keep bool) {
		if cur != tre.Root() {
			depths[cur.Id()] = depths[prev.Id()] + 1
		}
		return true
	})
	return depths
}

// Count leaves below each node
func countLeavesBelow(tre *tree.Tree, children [][]*tree.Node) []uint64 {
	below := make([]uint64, len(tre.Nodes()))
	tre.PostOrder(func(cur, prev *tree.Node, e *tree.Edge) (keep bool) {
		if cur.Tip() {
			below[cur.Id()] = 1
		} else {
			for _, c := range children[cur.Id()] {
				below[cur.Id()] += below[c.Id()]
			}
		}
		return true
	})
	return below
}

// Maps quartets to vertices where at least 3 taxa from the quartet exist below the vertex
func mapQuartetsToVertices(tre *tree.Tree, qCounts map[Quartet]uint32, leafsets []*bitset.BitSet) [][]Quartet {
	qSets := make([][]Quartet, len(tre.Nodes()))
	n, err := tre.NbTips()
	if err != nil {
		panic(err)
	}
	tre.PostOrder(func(cur, prev *tree.Node, e *tree.Edge) (keep bool) {
		qSets[cur.Id()] = make([]Quartet, 0)
		for q := range qCounts {
			found := 0
			for i := range 4 {
				if q.Taxon(i) >= uint16(n) {
					panic("cannot map quartet taxa to constraint tree")
				} else if leafsets[cur.Id()].Test(uint(q.Taxon(i))) {
					// } else if leafsets[cur.Id()][q.Taxa[i]] {
					found++
				}
			}
			if found >= 3 {
				qSets[cur.Id()] = append(qSets[cur.Id()], q)
			}
		}
		return true
	})
	return qSets
}

func makeTipIndexMap(tre *tree.Tree) map[uint16]int {
	tips := tre.Tips()
	tipMap := make(map[uint16]int, len(tips))
	for _, t := range tips {
		tipMap[uint16(t.TipIndex())] = t.Id()
	}
	return tipMap
}

// n2 (id) is in the leafset of n1 (id)
func (td *TreeData) InLeafset(n1ID, n2ID uint16) bool {
	return td.leafsets[n1ID].Test(uint(n2ID))
}

// Takes in the node ids of two nodes and returns the id of the LCA
func (td *TreeData) LCA(n1ID, n2ID int) int {
	return td.lca[n1ID][n2ID]
}

// Finds node's sibling -- assumes binary tree
func (td *TreeData) Sibling(node *tree.Node) *tree.Node {
	// assumes binary tree
	p, err := node.Parent()
	if err != nil {
		panic(err)
	}
	for _, c := range td.Children[p.Id()] {
		if c != node {
			return c
		}
	}
	panic("failed to find node sibling")
}

// Returns leafset as string for printing/testing
func (td *TreeData) LeafsetAsString(n *tree.Node) string {
	result := "{"
	tips := td.AllTipNames()
	for i := range len(tips) {
		// for i, t := range td.leafsets[n.Id()] {
		if td.leafsets[n.Id()].Test(uint(i)) {
			result += tips[i] + ","
		}
	}
	return result[:len(result)-1] + "}"
}

func (td *TreeData) TipToNodeID(idx uint16) int {
	return td.tipIndexMap[idx]
}

// Get quartets corresponding to a given node (by id)
func (td *TreeData) Quartets(nid int) []Quartet {
	if td.quartetSet == nil {
		panic("quartet set never initialized")
	}
	return td.quartetSet[nid]
}

// Get count of quartets with a particular topology
func (td *TreeData) NumQuartet(q Quartet) uint32 {
	if td.quartetSet == nil {
		panic("quartet counts never initialized")
	}
	return (*td.quartetCounts)[q]
}

// n2 is under n1
func (td *TreeData) Under(n1ID, n2ID int) bool {
	return td.LCA(n1ID, n2ID) == n1ID && n1ID != n2ID
}

// returns total number of quartets (all topologies)
func (td *TreeData) TotalNumQuartets() uint32 {
	var result uint32
	for _, count := range *td.quartetCounts {
		result += count
	}
	return result
}

func (td *TreeData) Clone() *TreeData {
	tre := td.Tree.Clone()
	return &TreeData{
		Tree:        *tre,
		Children:    children(tre),
		IdToNodes:   mapIdToNodes(tre),
		Depths:      td.Depths,
		leafsets:    td.leafsets,
		lca:         td.lca,
		tipIndexMap: td.tipIndexMap,
		NLeaves:     td.NLeaves,
	}
}
