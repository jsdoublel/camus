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

type Branch struct {
	IDs [2]int
}

func (br *Branch) Empty() bool {
	return br.IDs == [2]int{0, 0}
}

// Makes extended newick network out of newick tree and branch data computed by
// the CAMUS algorithm
func MakeNetwork(td *TreeData, branches []Branch) *Network {
	ret := make(map[string][2]int)
	for i, branch := range branches {
		ret[fmt.Sprintf("#H%d", i)] = branch.IDs // TODO: finish branch refactor
		u, w := td.IdToNodes[branch.IDs[Ui]], td.IdToNodes[branch.IDs[Wi]]
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

func (ntw *Network) Level1(td *TreeData) bool {
	branchs := make([]string, 0)
	for k := range ntw.Reticulations {
		branchs = append(branchs, k)
	}
	for i := range branchs {
		for j := i + 1; j < len(branchs); j++ {
			r1 := ntw.Reticulations[branchs[i]]
			r2 := ntw.Reticulations[branchs[j]]
			vR1 := td.LCA(r1[0], r1[1])
			vR2 := td.LCA(r2[0], r2[1])
			if vR1 == vR2 || illSorted(vR1, vR2, r1, td) || illSorted(vR2, vR1, r2, td) {
				return false
			}
		}
	}
	return true
}

func illSorted(v1, v2 int, r1 [2]int, td *TreeData) bool {
	return td.Under(v1, v2) && (td.Under(v2, r1[0]) || td.Under(v2, r1[1]))
}
