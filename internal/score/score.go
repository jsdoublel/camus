// Package implementing scoring for networks
package score

import (
	"errors"
	"fmt"
	"math"

	"github.com/evolbioinfo/gotree/tree"

	gr "github.com/jsdoublel/camus/internal/graphs"
	pr "github.com/jsdoublel/camus/internal/prep"
)

var ErrNotLevel1 = errors.New("not level-1")

// nodes needed for scoring reticulation
type reticulation struct {
	u    *tree.Node
	w    *tree.Node
	v    *tree.Node
	wSub *tree.Node
}

func ReticulationScore(ntw *gr.Network, gtrees []*tree.Tree) ([]*map[string]float64, error) {
	td := gr.MakeTreeData(ntw.NetTree, nil)
	if !ntw.Level1(td) {
		return nil, fmt.Errorf("network is %w", ErrNotLevel1)
	}
	reticulations := *getReticulationNodes(ntw, td)
	results := make([]*map[string]float64, len(gtrees))
	for i, gtre := range gtrees {
		if err := gtre.UpdateTipIndex(); err != nil {
			return nil, fmt.Errorf("gene tree %w", pr.ErrMulTree)
		}
		totals := make(map[string]uint)
		supported := make(map[string]uint)
		gtre.UnRoot()
		constMap, err := gr.MapIDsFromConstTree(gtre, ntw.NetTree)
		if err != nil {
			return nil, err
		}
		gtre.Quartets(false, func(q *tree.Quartet) {
			for label, branch := range reticulations {
				comp := QuartetScore(
					gr.QuartetFromTreeQ(q, constMap),
					branch.u,
					branch.w,
					branch.v,
					branch.wSub,
					td,
				)
				if comp != gr.Qdiff {
					totals[label] += 1
				}
				if comp == gr.Qeq {
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
func getReticulationNodes(ntw *gr.Network, td *gr.TreeData) *map[string]reticulation {
	result := make(map[string]reticulation)
	for label, branch := range ntw.Reticulations {
		uId, wId := branch.IDs[gr.Ui], branch.IDs[gr.Wi]
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
