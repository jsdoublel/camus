package prep

import (
	"fmt"

	"github.com/evolbioinfo/gotree/tree"
)

type TreeData struct {
	Tree        *tree.Tree     // Tree object
	Root        *tree.Node     // Root for which data is calculated
	Children    [][]*tree.Node // Children for each node
	Leafsets    [][]bool       // Leaves under each node
	IdToNodes   []*tree.Node   // Mapping between id and node pointer
	QuartetSet  [][]*Quartet   // Quartets relevant for each subtree
	leafRep     []*tree.Node   // A single leaf that can be looked up for each node (used for LCAs between internal nodes)
	lca         [][]int        // LCA for each pair of node id
	tipIndexMap map[int]int
}

func PreprocessTreeData(tre *tree.Tree, quartets []*Quartet) *TreeData {
	root := tre.Root()
	children := children(tre)
	lca, leafsets, leafReps := lcaAndLeafset(tre, children)
	idMap := mapIdToNodes(tre)
	quartetSets := mapQuartetsToVertices(tre, quartets, leafsets)
	tipIndexMap := makeTipIndexMap(tre)
	return &TreeData{Tree: tre, Root: root, Children: children, lca: lca,
		Leafsets: leafsets, leafRep: leafReps, IdToNodes: idMap,
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
	// return root == td.Root
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

/* panics (array out of bounds) if tree is not binary */
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

/* calculates the LCA for all pairs of leaves as well as the leaf set for every node */
func lcaAndLeafset(tre *tree.Tree, children [][]*tree.Node) ([][]int, [][]bool, []*tree.Node) {
	nLeaves := len(tre.Tips())
	nNodes := len(tre.Nodes())
	lca := make([][]int, nLeaves)
	for i := range nLeaves {
		lca[i] = make([]int, nLeaves)
	}
	leafset := make([][]bool, nNodes)
	repLeaves := make([]*tree.Node, nNodes)
	tre.PostOrder(func(cur, prev *tree.Node, e *tree.Edge) (keep bool) {
		leafset[cur.Id()] = make([]bool, nLeaves)
		if cur.Tip() {
			leafset[cur.Id()][cur.TipIndex()] = true
			lca[cur.TipIndex()][cur.TipIndex()] = cur.Id()
			repLeaves[cur.Id()] = cur
		} else {
			for i := 0; i < nLeaves; i++ {
				leafset[cur.Id()][i] = leafset[children[cur.Id()][0].Id()][i] || leafset[children[cur.Id()][1].Id()][i]
				for j := 0; j < nLeaves; j++ {
					if leafset[children[cur.Id()][0].Id()][i] == leafset[children[cur.Id()][1].Id()][j] {
						lca[i][j] = cur.Id()
						lca[j][i] = cur.Id()
					}
				}
			}
			repLeaves[cur.Id()] = repLeaves[children[cur.Id()][0].Id()]
		}
		return true
	})
	return lca, leafset, repLeaves
}

func mapQuartetsToVertices(tre *tree.Tree, quartets []*Quartet, leafsets [][]bool) [][]*Quartet {
	quartetSets := make([][]*Quartet, len(tre.Nodes()))
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
		tipMap[t.Id()] = t.TipIndex()
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

/* takes in the node ids of two nodes (not the tip indices) and returns the id of the LCA */
func (td *TreeData) LCA(n1ID, n2ID int) int {
	n1, n2 := td.IdToNodes[n1ID], td.IdToNodes[n2ID]
	if !n1.Tip() {
		n1 = td.leafRep[n1ID]
	}
	if !n2.Tip() {
		n2 = td.leafRep[n2ID]
	}
	fmt.Printf("lca %s %s\n", n1.Name(), n2.Name())
	return td.lca[td.tipIndexMap[n1.Id()]][td.tipIndexMap[n2.Id()]]
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
