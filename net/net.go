package net

import (
	"fmt"
	"strings"

	"github.com/evolbioinfo/gotree/tree"

	"github.com/jsdoublel/camus/prep"
)

type Network struct {
	ExTree        *tree.Tree // tree from extended newick
	Reticulations [][2]int   // reticulation branches
}

// Makes extended newick network out of newick tree and branch data in somewhat
// hacky way.
func MakeNetwork(td *prep.TreeData, branches [][2]int) *Network {
	for i, branch := range branches {
		u, w := td.IdToNodes[branch[0]], td.IdToNodes[branch[1]]
		uEdge, err := u.ParentEdge()
		if err != nil {
			panic(fmt.Sprintf("error in MakeNetwork getting u (id %d): %s", u.Id(), err))
		}
		r := td.Tree.NewNode()
		r.SetName(fmt.Sprintf("#H%d", i))
		td.Tree.GraftTipOnEdge(r, uEdge)
		r = td.Tree.NewNode()
		r.SetName("####")
		wEdge, err := w.ParentEdge()
		if err != nil {
			panic(fmt.Sprintf("error in MakeNetwork getting w: %s", err))
		}
		td.Tree.GraftTipOnEdge(r, wEdge)
		p, err := r.Parent()
		if err != nil {
			panic(fmt.Sprintf("error in MakeNetwork after grafting w: %s", err))
		}
		p.SetName(fmt.Sprintf("#H%d", i))
	}
	cleanTree(td.Tree)
	return &Network{ExTree: td.Tree, Reticulations: branches}
}

func (ntw *Network) Newick() string {
	nwk := ntw.ExTree.Newick()
	nwk = strings.ReplaceAll(nwk, "####,", "")
	nwk = strings.ReplaceAll(nwk, ",####", "")
	return nwk
}

// Deletes all branch lengths and support values (since they might be missleading)
func cleanTree(tre *tree.Tree) {
	tre.PostOrder(func(cur, prev *tree.Node, e *tree.Edge) (keep bool) {
		if e != nil {
			e.SetSupport(tree.NIL_SUPPORT)
			e.SetLength(tree.NIL_LENGTH)
		}
		return true
	})
}
