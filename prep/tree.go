package prep

import (
	"fmt"

	"github.com/evolbioinfo/gotree/tree"
)

/* panics (array out of bounds) if tree is not binary */
func GetChildren(node *tree.Node) []*tree.Node {
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

func IsBinary(tree *tree.Tree) bool {
	root := tree.Root()
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
	children := GetChildren(node)
	if len(children) != 2 {
		return false
	}
	return isBinaryRecusive(children[0]) && isBinaryRecusive(children[1])
}
