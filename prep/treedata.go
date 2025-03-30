package prep

import (
	"github.com/bits-and-blooms/bitset"
	"github.com/evolbioinfo/gotree/tree"

	"github.com/jsdoublel/camus/qrt"
)

type TreeData struct {
	Tree          *tree.Tree            // Tree object
	Root          *tree.Node            // Root for which data is calculated
	Children      [][]*tree.Node        // Children for each node
	IdToNodes     []*tree.Node          // Mapping between id and node pointer
	QuartetSet    [][]*qrt.Quartet      // Quartets relevant for each subtree
	QuartetCounts *map[qrt.Quartet]uint // Count of each unqiue quartet topology
	Depths        []int                 // Distance from all nodes to the root
	leafsets      []*bitset.BitSet      // Leaves under each node
	lca           [][]int               // LCA for each pair of node id
	tipIndexMap   map[int]int           // Tip index to node id map
	NLeaves       int                   // Number of leaves
}

func MakeTreeData(tre *tree.Tree, qCounts *map[qrt.Quartet]uint) *TreeData {
	root := tre.Root()
	children := children(tre)
	leafsets := calcLeafset(tre, children)
	lca := calcLCAs(tre, children)
	depths := calcDepths(tre)
	idMap := mapIdToNodes(tre)
	qSets := mapQuartetsToVertices(tre, qCounts, leafsets)
	tipIndexMap := makeTipIndexMap(tre)
	return &TreeData{Tree: tre,
		Root:          root,
		Children:      children,
		lca:           lca,
		leafsets:      leafsets,
		IdToNodes:     idMap,
		Depths:        depths,
		QuartetSet:    qSets,
		QuartetCounts: qCounts,
		tipIndexMap:   tipIndexMap,
		NLeaves:       len(tre.AllTipNames()),
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
	root := td.Tree.Root()
	if root != td.Root {
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
			children[cur.Id()] = getChildren(cur)
		}
		return true
	})
	return children
}

// Get children of node
func getChildren(node *tree.Node) []*tree.Node {
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
			leafset[cur.Id()] = leafset[children[cur.Id()][0].Id()].Union(leafset[children[cur.Id()][1].Id()])
		}
		return true
	})
	return leafset
}

// Caculates the LCA for every pair of nodes
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
			leftId, rightId := children[cur.Id()][0].Id(), children[cur.Id()][1].Id()
			for i := range nNodes {
				below[cur.Id()][i] = below[leftId][i] || below[rightId][i] || i == cur.Id()
			}
			for i := range nNodes {
				for j := range nNodes {
					if below[leftId][i] && below[rightId][j] ||
						i == cur.Id() && below[rightId][j] ||
						i == cur.Id() && below[leftId][j] {
						lca[i][j] = cur.Id()
						lca[j][i] = cur.Id()
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

// Maps quartets to vertices where at least 3 taxa from the quartet exist below the vertex
func mapQuartetsToVertices(tre *tree.Tree, qCounts *map[qrt.Quartet]uint, leafsets []*bitset.BitSet) [][]*qrt.Quartet {
	qSets := make([][]*qrt.Quartet, len(tre.Nodes()))
	n, err := tre.NbTips()
	if err != nil {
		panic(err)
	}
	tre.PostOrder(func(cur, prev *tree.Node, e *tree.Edge) (keep bool) {
		qSets[cur.Id()] = make([]*qrt.Quartet, 0)
		for q := range *qCounts {
			found := 0
			for i := range 4 {
				if q.Taxa[i] >= n {
					panic("cannot map quartet taxa to constraint tree")
				} else if leafsets[cur.Id()].Test(uint(q.Taxa[i])) {
					// } else if leafsets[cur.Id()][q.Taxa[i]] {
					found++
				}
			}
			if found >= 3 {
				qSets[cur.Id()] = append(qSets[cur.Id()], &q)
			}
		}
		return true
	})
	return qSets
}

func makeTipIndexMap(tre *tree.Tree) map[int]int {
	tips := tre.Tips()
	tipMap := make(map[int]int, len(tips))
	for _, t := range tips {
		tipMap[t.TipIndex()] = t.Id()
	}
	return tipMap
}

func (td *TreeData) InLeafset(n1ID, n2ID int) bool {
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
	tips := td.Tree.AllTipNames()
	for i := range len(tips) {
		// for i, t := range td.leafsets[n.Id()] {
		if td.leafsets[n.Id()].Test(uint(i)) {
			result += tips[i] + ","
		}
	}
	return result[:len(result)-1] + "}"
}

func (td *TreeData) NodeID(idx int) int {
	return td.tipIndexMap[idx]
}
