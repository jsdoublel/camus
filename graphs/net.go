package graphs

import (
	"fmt"
	"slices"
	"strings"

	"github.com/evolbioinfo/gotree/tree"
)

type Network struct {
	NetTree       *tree.Tree        // tree from extended newick
	Reticulations map[string]Branch // reticulation branches
}

const (
	Ui = iota // index of u and w in reticulation array
	Wi
)

type Branch struct {
	IDs [2]int // {0: u, 1: w}
}

func (br Branch) Empty() bool {
	return br.IDs == [2]int{0, 0}
}

func (br Branch) Collide(br2 Branch) bool {
	return (br.IDs[0] == br2.IDs[0] ||
		br.IDs[0] == br2.IDs[1] ||
		br.IDs[1] == br2.IDs[0] ||
		br.IDs[1] == br2.IDs[1])
}

// Makes extended newick network out of newick tree and branch data computed by
// the CAMUS algorithm
func MakeNetwork(td *TreeData, branches []Branch) *Network {
	td = td.Clone()
	ret := make(map[string]Branch)
	slices.SortFunc(branches, func(br1, br2 Branch) int {
		if br1.Collide(br2) {
			if td.Under(br1.IDs[0], br2.IDs[0]) ||
				td.Under(br1.IDs[0], br2.IDs[1]) ||
				td.Under(br1.IDs[1], br2.IDs[0]) ||
				td.Under(br1.IDs[1], br2.IDs[1]) {
				return -1
			} else {
				return 1
			}
		}
		return 0
	})
	for i, branch := range branches {
		ret[fmt.Sprintf("#H%d", i)] = branch
		u, w := td.IdToNodes[branch.IDs[Ui]], td.IdToNodes[branch.IDs[Wi]]
		uEdge, err := u.ParentEdge()
		if err != nil {
			panic(fmt.Sprintf("error in MakeNetwork getting u (id %d): %s", u.Id(), err))
		}
		r := td.NewNode()
		r.SetName(fmt.Sprintf("#H%d", i))
		if _, _, _, err := td.GraftTipOnEdge(r, uEdge); err != nil {
			panic(err)
		}
		r = td.NewNode()
		r.SetName("####")
		wEdge, err := w.ParentEdge()
		if err != nil {
			panic(fmt.Sprintf("error in MakeNetwork getting w: %s", err))
		}
		if _, _, _, err := td.GraftTipOnEdge(r, wEdge); err != nil {
			panic(err)
		}
		p, err := r.Parent()
		if err != nil {
			panic(fmt.Sprintf("error in MakeNetwork after grafting w: %s", err))
		}
		p.SetName(fmt.Sprintf("#H%d", i))
	}
	cleanTree(&td.Tree)
	return &Network{NetTree: &td.Tree, Reticulations: ret}
}

func (ntw *Network) Newick() string {
	nwk := ntw.NetTree.Newick()
	nwk = strings.ReplaceAll(nwk, "####,", "")
	nwk = strings.ReplaceAll(nwk, ",####", "")
	return nwk
}

// Deletes all branch lengths and support values (since they might be misleading)
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
	branches := make([]string, 0)
	for k := range ntw.Reticulations {
		branches = append(branches, k)
	}
	for i := range branches {
		for j := i + 1; j < len(branches); j++ {
			r1 := ntw.Reticulations[branches[i]]
			r2 := ntw.Reticulations[branches[j]]
			vR1 := td.LCA(r1.IDs[0], r1.IDs[1])
			vR2 := td.LCA(r2.IDs[0], r2.IDs[1])
			if vR1 == vR2 || illSorted(vR1, vR2, r1, td) || illSorted(vR2, vR1, r2, td) {
				return false
			}
		}
	}
	return true
}

func illSorted(v1, v2 int, r1 Branch, td *TreeData) bool {
	return td.Under(v1, v2) && (td.Under(v2, r1.IDs[0]) || td.Under(v2, r1.IDs[1]))
}
