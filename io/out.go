package io

import (
	"camus/prep"
	"fmt"
	"strings"
)

func MakeNetwork(td *prep.TreeData, branches [][2]int) string {
	for i, branch := range branches {
		u, w := td.IdToNodes[branch[0]], td.IdToNodes[branch[1]]
		uEdge, err := u.ParentEdge()
		if err != nil {
			panic(err)
		}
		r := td.Tree.NewNode()
		r.SetName(fmt.Sprintf("#%d", i))
		td.Tree.GraftTipOnEdge(r, uEdge)
		r = td.Tree.NewNode()
		r.SetName("####")
		wEdge, err := w.ParentEdge()
		if err != nil {
			panic(err)
		}
		td.Tree.GraftTipOnEdge(r, wEdge)
		p, err := r.Parent()
		if err != nil {
			panic(err)
		}
		p.SetName(fmt.Sprintf("#%d", i))
	}
	return fixNetwork(td.Tree.Newick())
}

func fixNetwork(nwk string) string {
	nwk = strings.ReplaceAll(nwk, "####:1,", "")
	nwk = strings.ReplaceAll(nwk, ",####", "")
	return nwk
}
