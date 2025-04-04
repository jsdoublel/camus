package graphs

import (
	"fmt"
	"strings"

	"github.com/evolbioinfo/gotree/tree"
)

type Network struct {
	NetTree       *tree.Tree        // tree from extended newick
	Reticulations map[string][2]int // reticulation branches
}

const (
	Ui = iota // index of u and w in reticulation array
	Wi
)

// Makes extended newick network out of newick tree and branch data computed by
// the CAMUS algorithm
func MakeNetwork(td *TreeData, branches [][2]int) *Network {
	ret := make(map[string][2]int)
	for i, branch := range branches {
		ret[fmt.Sprintf("#H%d", i)] = branch
		u, w := td.IdToNodes[branch[Ui]], td.IdToNodes[branch[Wi]]
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
	return &Network{NetTree: td.Tree, Reticulations: ret}
}

func (ntw *Network) Newick() string {
	nwk := ntw.NetTree.Newick()
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
