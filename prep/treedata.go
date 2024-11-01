package prep

import (
	"fmt"

	"github.com/evolbioinfo/gotree/tree"
)

type TreeData struct {
	Tree        *tree.Tree     // Tree object
	Root        *tree.Node     // Root for which data is calculated
	Children    [][]*tree.Node // Children for each node
	IdToNodes   []*tree.Node   // Mapping between id and node pointer
	QuartetSet  [][]*Quartet   // Quartets relevant for each subtree
	Depths      []int          // Distance from all nodes to the root (which as it turns out is the correct definition gotree!!!)
	leafsets    [][]bool       // Leaves under each node
	lca         [][]int        // LCA for each pair of node id
	tipIndexMap map[int]int
}

func PreprocessTreeData(tre *tree.Tree, quartets []*Quartet) *TreeData {
	root := tre.Root()
	children := children(tre)
	// lca, leafsets, leafReps := lcaAndLeafset(tre, children)
	leafsets := calcLeafset(tre, children)
	lca := calcLCAs(tre, children)
	depths := calcDepths(tre)
	idMap := mapIdToNodes(tre)
	quartetSets := mapQuartetsToVertices(tre, quartets, leafsets)
	tipIndexMap := makeTipIndexMap(tre)
	return &TreeData{Tree: tre, Root: root, Children: children, lca: lca,
		leafsets: leafsets, IdToNodes: idMap, Depths: depths,
		QuartetSet: quartetSets, tipIndexMap: tipIndexMap}
}

func mapIdToNodes(tre *tree.Tree) []*tree.Node {
	idMap := make([]*tree.Node, len(tre.Nodes()))
	tre.PostOrder(func(cur, prev *tree.Node, e *tree.Edge) (keep bool) {
		idMap[cur.Id()] = cur
		return true
	})
	return idMap
}

/* verify that tree still has the same root, and thus the data is still applicable */
func (td *TreeData) Verify() {
	root := td.Tree.Root()
	if root != td.Root {
		panic("TreeData root is wrong!")
	}
}

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

/* can panic (array out of bounds) if tree is not binary */
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

/* calculates the leafset for every node */
func calcLeafset(tre *tree.Tree, children [][]*tree.Node) [][]bool {
	nLeaves, nNodes := len(tre.Tips()), len(tre.Nodes())
	leafset := make([][]bool, nNodes)
	tre.PostOrder(func(cur, prev *tree.Node, e *tree.Edge) (keep bool) {
		leafset[cur.Id()] = make([]bool, nLeaves)
		if cur.Tip() {
			leafset[cur.Id()][cur.TipIndex()] = true
		} else {
			for i := 0; i < nLeaves; i++ {
				leafset[cur.Id()][i] = leafset[children[cur.Id()][0].Id()][i] || leafset[children[cur.Id()][1].Id()][i]
			}
		}
		return true
	})
	return leafset
}

/* caculates the LCA for every pair of nodes */
func calcLCAs(tre *tree.Tree, children [][]*tree.Node) [][]int {
	nNodes := len(tre.Nodes())
	lca, below := make([][]int, nNodes), make([][]bool, nNodes) // below[i][j] = true means node j is below node i
	for i := range nNodes {
		lca[i] = make([]int, nNodes)
		below[i] = make([]bool, nNodes)
	}
	tre.PostOrder(func(cur, prev *tree.Node, e *tree.Edge) (keep bool) {
		below[cur.Id()][cur.Id()] = true
		if cur.Tip() {
			lca[cur.Id()][cur.Id()] = cur.Id()
		} else {
			for i := range nNodes {
				below[cur.Id()][i] = below[children[cur.Id()][0].Id()][i] || below[children[cur.Id()][1].Id()][i]
				for j := range nNodes {
					if below[children[cur.Id()][0].Id()][i] && below[children[cur.Id()][1].Id()][j] {
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

// /* calculates the LCA for all pairs of leaves as well as the leaf set for every node */
// func lcaAndLeafset(tre *tree.Tree, children [][]*tree.Node) ([][]int, [][]bool, []*tree.Node) {
// 	nLeaves := len(tre.Tips())
// 	nNodes := len(tre.Nodes())
// 	lca := make([][]int, nLeaves)
// 	for i := range nLeaves {
// 		lca[i] = make([]int, nLeaves)
// 	}
// 	leafset := make([][]bool, nNodes)
// 	repLeaves := make([]*tree.Node, nNodes)
// 	tre.PostOrder(func(cur, prev *tree.Node, e *tree.Edge) (keep bool) {
// 		leafset[cur.Id()] = make([]bool, nLeaves)
// 		if cur.Tip() {
// 			leafset[cur.Id()][cur.TipIndex()] = true
// 			lca[cur.TipIndex()][cur.TipIndex()] = cur.Id()
// 			repLeaves[cur.Id()] = cur
// 		} else {
// 			for i := 0; i < nLeaves; i++ {
// 				leafset[cur.Id()][i] = leafset[children[cur.Id()][0].Id()][i] || leafset[children[cur.Id()][1].Id()][i]
// 				for j := 0; j < nLeaves; j++ {
// 					if leafset[children[cur.Id()][0].Id()][i] && leafset[children[cur.Id()][1].Id()][j] {
// 						lca[i][j] = cur.Id()
// 						lca[j][i] = cur.Id()
// 					}
// 				}
// 			}
// 			repLeaves[cur.Id()] = repLeaves[children[cur.Id()][0].Id()]
// 		}
// 		return true
// 	})
// 	return lca, leafset, repLeaves
// }

func mapQuartetsToVertices(tre *tree.Tree, quartets []*Quartet, leafsets [][]bool) [][]*Quartet {
	quartetSets := make([][]*Quartet, len(tre.Nodes()))
	// for i := range len(tre.Nodes()) {
	// 	quartetSets[i] = quartets
	// }
	n := len(tre.Tips())
	tre.PostOrder(func(cur, prev *tree.Node, e *tree.Edge) (keep bool) {
		quartetSets[cur.Id()] = make([]*Quartet, 0)
		for _, q := range quartets {
			add := true
			for i := 0; i < 4; i++ {
				if q.Taxa[i] >= n {
					panic("cannot map quartet taxa to constraint tree")
				} else if !leafsets[cur.Id()][q.Taxa[i]] {
					add = false
					break
				}
			}
			if add {
				quartetSets[cur.Id()] = append(quartetSets[cur.Id()], q)
			}
		}
		return true
	})
	return quartetSets
}

func makeTipIndexMap(tre *tree.Tree) map[int]int {
	tips := tre.Tips()
	tipMap := make(map[int]int, len(tips))
	for _, t := range tips {
		// tipMap[t.Id()] = t.TipIndex()
		tipMap[t.TipIndex()] = t.Id()
	}
	return tipMap
}

func IsBinary(tre *tree.Tree) bool {
	root := tre.Root()
	neighbors := root.Neigh()
	if len(neighbors) != 2 {
		return false
	}
	return isBinaryRecusive(neighbors[0]) && isBinaryRecusive(neighbors[1])
}

func isBinaryRecusive(node *tree.Node) bool {
	if node.Tip() {
		return true
	}
	children := getChildren(node)
	if len(children) != 2 {
		return false
	}
	return isBinaryRecusive(children[0]) && isBinaryRecusive(children[1])
}

func (td *TreeData) InLeafset(n1ID, n2ID int) bool {
	return td.leafsets[n1ID][n2ID]
}

/* takes in the node ids of two nodes and returns the id of the LCA */
func (td *TreeData) LCA(n1ID, n2ID int) int {
	// n1 := td.IdToNodes[n1ID]
	// n1IDtip := n1.TipIndex()
	// if !n1.Tip() {
	// 	n1IDtip = td.leafRep[n1ID].TipIndex()
	// }
	return td.lca[n1ID][n2ID]
}

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
