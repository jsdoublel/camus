package prep

import (
	"fmt"

	"github.com/evolbioinfo/gotree/tree"
)

type TreeData struct {
	Tree        *tree.Tree
	Root        *tree.Node
	Children    [][]*tree.Node
	LCA         [][]uint
	Leafsets    [][]bool
	QuartetSets [][]*Quartet
}

func PreprocessTreeData(tre *tree.Tree, quartets []*Quartet) *TreeData {
	root := tre.Root()
	children := children(tre)
	lca, leafsets := lcaAndLeafset(tre, children)
	quartetSets := mapQuartetsToVertices(tre, quartets, leafsets)
	return &TreeData{Tree: tre, Root: root, Children: children, LCA: lca, Leafsets: leafsets, QuartetSets: quartetSets}
}

/* verify that tree still has the same root, and thus the data is still applicable */
func (td *TreeData) Verify() bool {
	root := td.Tree.Root()
	return root == td.Root
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
func lcaAndLeafset(tre *tree.Tree, children [][]*tree.Node) ([][]uint, [][]bool) {
	nLeaves := len(tre.Tips())
	nNodes := len(tre.Nodes())
	lca := make([][]uint, nLeaves)
	for i := range nLeaves {
		lca[i] = make([]uint, nLeaves)
	}
	leafset := make([][]bool, nNodes)
	tre.PostOrder(func(cur, prev *tree.Node, e *tree.Edge) (keep bool) {
		leafset[cur.Id()] = make([]bool, nLeaves)
		if cur.Tip() {
			leafset[cur.Id()][cur.TipIndex()] = true
			lca[cur.TipIndex()][cur.TipIndex()] = uint(cur.Id())
		} else {
			for i := 0; i < nLeaves; i++ {
				leafset[cur.Id()][i] = leafset[children[cur.Id()][0].Id()][i] || leafset[children[cur.Id()][1].Id()][i]
				for j := 0; j < nLeaves; j++ {
					if leafset[children[cur.Id()][0].Id()][i] == leafset[children[cur.Id()][1].Id()][j] {
						lca[i][j] = uint(cur.Id())
						lca[j][i] = uint(cur.Id())
					}
				}
			}
		}
		return true
	})
	return lca, leafset
}

func mapQuartetsToVertices(tre *tree.Tree, quartets []*Quartet, leafsets [][]bool) [][]*Quartet {
	quartetSets := make([][]*Quartet, len(tre.Nodes()))
	n := len(tre.Tips())
	tre.PostOrder(func(cur, prev *tree.Node, e *tree.Edge) (keep bool) {
		quartetSets[cur.Id()] = make([]*Quartet, 0)
		for _, q := range quartets {
			add := true
			for i := 0; i < 4; i++ {
				if q.taxa[i] >= uint(n) {
					panic("cannot map quartet taxa to constraint tree")
				} else if !leafsets[cur.Id()][q.taxa[i]] {
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
