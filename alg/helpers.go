package alg

import "github.com/evolbioinfo/gotree/tree"

func SubtreePostOrder(cur *tree.Node, f func(cur, otherSubtree *tree.Node)) {
	if !cur.Tip() {
		children := make([]*tree.Node, 0)
		for _, n := range cur.Neigh() {
			if p, err := cur.Parent(); err != nil {
				panic(err)
			} else if n != p {
				children = append(children, n)
			}
		}
		if len(children) != 2 {
			panic("tree is not binary")
		}
		subtreePostOrderHelper(children[0], children[1], f)
		subtreePostOrderHelper(children[1], children[0], f)
	}
}

func subtreePostOrderHelper(cur, otherSubtree *tree.Node, f func(cur, otherSubtree *tree.Node)) {
	for _, n := range cur.Neigh() {
		if p, err := cur.Parent(); err != nil {
			panic(err)
		} else if n != p {
			subtreePostOrderHelper(n, otherSubtree, f)
		}
	}
	f(cur, otherSubtree)
}

func SubtreePreOrder(cur *tree.Node, f func(cur *tree.Node)) {
	f(cur)
	for _, n := range cur.Neigh() {
		if p, err := cur.Parent(); err != nil {
			panic(err)
		} else if n != p {
			SubtreePreOrder(n, f)
		}
	}
}
