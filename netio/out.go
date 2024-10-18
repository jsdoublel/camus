package netio

import (
	"camus/prep"
	"fmt"
	"strings"

	"github.com/evolbioinfo/gotree/tree"
)

func MakeNetwork(td *prep.TreeData, branches [][2]int) string {
	for i, branch := range branches {
		u, w := td.IdToNodes[branch[0]], td.IdToNodes[branch[1]]
		uEdge, err := u.ParentEdge()
		if err != nil {
			panic(fmt.Sprintf("error in MakeNetwork getting u (id %d): %s", u.Id(), err))
		}
		r := td.Tree.NewNode()
		r.SetName(fmt.Sprintf("#%d", i))
		td.Tree.GraftTipOnEdge(r, uEdge)
		r = td.Tree.NewNode()
		r.SetName("####")
		wEdge, err := w.ParentEdge()
		if err != nil {
			panic(fmt.Sprintf("error in MakeNetwork getting w: %s", err))
		}
		td.Tree.GraftTipOnEdge(r, wEdge)
		td.Tree.ReinitIndexes()
		p, err := r.Parent()
		if err != nil {
			panic(fmt.Sprintf("error in MakeNetwork after grafting w: %s", err))
		}
		p.SetName(fmt.Sprintf("#%d", i))
	}
	deleteAllBranchLengths(td.Tree)
	return fixNetwork(td.Tree.Newick())
}

func fixNetwork(nwk string) string {
	nwk = strings.ReplaceAll(nwk, "####,", "")
	nwk = strings.ReplaceAll(nwk, ",####", "")
	return nwk
}

func deleteAllBranchLengths(tre *tree.Tree) {
	tre.PostOrder(func(cur, prev *tree.Node, e *tree.Edge) (keep bool) {
		if e != nil {
			e.SetLength(tree.NIL_LENGTH)
		}
		return true
	})
}
