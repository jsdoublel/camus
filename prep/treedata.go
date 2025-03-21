package prep

import (
	"fmt"

	"github.com/evolbioinfo/gotree/tree"
)

type TreeData struct {
	Tree          *tree.Tree       // Tree object
	Root          *tree.Node       // Root for which data is calculated
	Children      [][]*tree.Node   // Children for each node
	IdToNodes     []*tree.Node     // Mapping between id and node pointer
	QuartetSet    [][]*Quartet     // Quartets relevant for each subtree
	QuartetCounts map[Quartet]uint // Count of each unqiue quartet topology
	Depths        []int            // Distance from all nodes to the root
	leafsets      [][]bool         // Leaves under each node
	lca           [][]int          // LCA for each pair of node id
	tipIndexMap   map[int]int
}

func MakeTreeData(tre *tree.Tree, qCounts map[Quartet]uint) *TreeData {
	root := tre.Root()
	children := children(tre)
	leafsets := calcLeafset(tre, children)
	lca := calcLCAs(tre, children)
	depths := calcDepths(tre)
	idMap := mapIdToNodes(tre)
	qSets := mapQuartetsToVertices(tre, qCounts, leafsets)
	tipIndexMap := makeTipIndexMap(tre)
	return &TreeData{Tree: tre, Root: root, Children: children, lca: lca,
		leafsets: leafsets, IdToNodes: idMap, Depths: depths,
		QuartetSet: qSets, QuartetCounts: qCounts, tipIndexMap: tipIndexMap}
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

// Get children of node -- can panic (array out of bounds) if tree is not binary
func getChildren(node *tree.Node) []*tree.Node {
	children := make([]*tree.Node, 2)
	p, err := node.Parent()
	if err != nil && err.Error() == "The node has more than one parent" {
		panic(fmt.Errorf("GetChildren %w", err))
	}
	i := 0
	for _, u := range node.Neigh() {
		if u != p {
			children[i] = u
			i++
		}
	}
	return children
}

// Calculates the leafset for every node
func calcLeafset(tre *tree.Tree, children [][]*tree.Node) [][]bool {
	nLeaves, nNodes := len(tre.Tips()), len(tre.Nodes())
	leafset := make([][]bool, nNodes)
	tre.PostOrder(func(cur, prev *tree.Node, e *tree.Edge) (keep bool) {
		leafset[cur.Id()] = make([]bool, nLeaves)
		if cur.Tip() {
			leafset[cur.Id()][cur.TipIndex()] = true
		} else {
			for i := range nLeaves {
				leafset[cur.Id()][i] = leafset[children[cur.Id()][0].Id()][i] || leafset[children[cur.Id()][1].Id()][i]
			}
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
func mapQuartetsToVertices(tre *tree.Tree, qCounts map[Quartet]uint, leafsets [][]bool) [][]*Quartet {
	qSets := make([][]*Quartet, len(tre.Nodes()))
	n := len(tre.Tips())
	tre.PostOrder(func(cur, prev *tree.Node, e *tree.Edge) (keep bool) {
		qSets[cur.Id()] = make([]*Quartet, 0)
		for q := range qCounts {
			found := 0
			for i := range 4 {
				if q.Taxa[i] >= n {
					panic("cannot map quartet taxa to constraint tree")
				} else if leafsets[cur.Id()][q.Taxa[i]] {
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
	return td.leafsets[n1ID][n2ID]
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

func (td *TreeData) NLeaves() int {
	return len(td.leafsets[0])
}

// Returns leafset as string for printing/testing
func (td *TreeData) LeafsetAsString(n *tree.Node) string {
	result := "{"
	tips := td.Tree.AllTipNames()
	for i, t := range td.leafsets[n.Id()] {
		if t {
			result += tips[i] + ","
		}
	}
	return result[:len(result)-1] + "}"
}

func (td *TreeData) NodeID(idx int) int {
	return td.tipIndexMap[idx]
}
