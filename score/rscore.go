package score

import (
	"fmt"
	"math"

	"github.com/evolbioinfo/gotree/tree"

	"github.com/jsdoublel/camus/graphs"
	"github.com/jsdoublel/camus/infer"
	"github.com/jsdoublel/camus/prep"
)

// nodes needed for scoring reticulation
type reticulation struct {
	u    *tree.Node
	w    *tree.Node
	v    *tree.Node
	wSub *tree.Node
}

func CalculateReticulationScore(ntw *graphs.Network, gtrees []*tree.Tree) ([]*map[string]float64, error) {
	td := graphs.MakeTreeData(ntw.NetTree, nil)
	reticulations := *getReticulationNodes(ntw, td)
	results := make([]*map[string]float64, len(gtrees))
	for i, gtre := range gtrees {
		prep.LogEveryNPercent(i, 10, len(gtrees), fmt.Sprintf("scoring gene tree %d of %d", i, len(gtrees)))
		gtre.UpdateTipIndex()
		totals := make(map[string]uint)
		supported := make(map[string]uint)
		gtre.UnRoot()
		constMap, err := graphs.MapIDsFromConstTree(gtre, ntw.NetTree)
		if err != nil {
			return nil, err
		}
		gtre.Quartets(false, func(q *tree.Quartet) {
			for label, branch := range reticulations {
				comp := infer.QuartetScore(
					graphs.QuartetFromTreeQ(q, constMap),
					branch.u,
					branch.w,
					branch.v,
					branch.wSub,
					td,
				)
				if comp != graphs.Qdiff {
					totals[label] += 1
				}
				if comp == graphs.Qeq {
					supported[label] += 1
				}
			}
		})
		gtreeResult := make(map[string]float64)
		for label := range reticulations {
			if totals[label] != 0 {
				gtreeResult[label] = float64(supported[label]) / float64(totals[label])
			} else {
				gtreeResult[label] = math.NaN()
			}
		}
		results[i] = &gtreeResult
	}
	return results, nil
}

// Get reticulation name to node map
func getReticulationNodes(ntw *graphs.Network, td *graphs.TreeData) *map[string]reticulation {
	result := make(map[string]reticulation)
	for label, branch := range ntw.Reticulations {
		uId, wId := branch[graphs.Ui], branch[graphs.Wi]
		vId := td.LCA(uId, wId)
		for _, neigh := range td.IdToNodes[vId].Neigh() {
			if td.LCA(vId, wId) == vId {
				result[label] = reticulation{
					u:    td.IdToNodes[uId],
					w:    td.IdToNodes[wId],
					v:    td.IdToNodes[vId],
					wSub: neigh,
				}
			}
		}
	}
	if len(result) != len(ntw.Reticulations) {
		panic(fmt.Sprintf("could not map reticulations to nodes %v", result))
	}
	return &result
}
