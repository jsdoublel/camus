package alg

import (
	"fmt"
	"os"

	"github.com/evolbioinfo/gotree/tree"

	"camus/prep"
)

type DP struct {
	DP        []uint         // score for each dp subproblem
	Branches  [][2]int       // branch for each dp subproblem
	Backtrace [][]int        // prev subproblems for each dp subproblem
	TreeData  *prep.TreeData // preprocessed data for our constraint tree
}

func CAMUS(tre *tree.Tree, geneTrees []*tree.Tree) ([][2]int, error) {
	fmt.Fprint(os.Stderr, "beginning data preprocessing")
	td, err := prep.Preprocess(tre, geneTrees)
	if err != nil {
		return nil, fmt.Errorf("Preprocess error: %w", err)
	}
	fmt.Fprint(os.Stderr, "preprocessing finished, beginning CAMUS")
	n := len(td.Tree.Nodes())
	dp := &DP{DP: make([]uint, n, n), Branches: make([][2]int, n), Backtrace: make([][]int, n), TreeData: td}
	return dp.RunDP(), nil
}

func (dp *DP) RunDP() [][2]int {
	dp.TreeData.Tree.PostOrder(func(cur, prev *tree.Node, e *tree.Edge) (keep bool) {
		if !cur.Tip() { // default value is 0, so we don't need to code a base case
			lID, rID := dp.TreeData.Children[cur.Id()][0].Id(), dp.TreeData.Children[cur.Id()][1].Id()
			score, branch := dp.score(cur)
			noEdgeScore := dp.DP[lID] + dp.DP[rID]
			if score > noEdgeScore {
				dp.DP[cur.Id()] = score
				dp.Branches[cur.Id()] = branch
			} else {
				dp.DP[cur.Id()] = noEdgeScore
				dp.Backtrace[cur.Id()] = []int{lID, rID}
			}
		}
		return true
	})
	result := dp.traceback()
	return result
}

/* calculates score for given top node v; returns score, best edge, and dp lookups*/
func (dp *DP) score(v *tree.Node) (uint, [2]int) {
	bestScore := uint(0)
	bestEdge := [2]int{0, 0}
	vPathDPScores := dp.accumlateDPScores(v)
	if v != dp.TreeData.Root {
		var wID int
		bestScore, wID = dp.scoreU(v, v, v, vPathDPScores)
		bestEdge = [2]int{v.Id(), wID}
	}
	SubtreePostOrder(v, func(cur, otherSubtree *tree.Node) {
		score, wID := dp.scoreU(cur, otherSubtree, v, vPathDPScores)
		if score >= bestScore {
			bestScore, bestEdge = score, [2]int{cur.Id(), wID}
		}
	})
	return bestScore, bestEdge
}

func (dp *DP) accumlateDPScores(v *tree.Node) map[int]uint {
	pathScores := make(map[int]uint)
	SubtreePreOrder(v, func(cur *tree.Node) {
		if p, err := cur.Parent(); v != dp.TreeData.Root && err != nil {
			panic(err)
		} else if p != v {
			pathScores[cur.Id()] = pathScores[p.Id()] + dp.DP[dp.TreeData.Sibling(cur).Id()]
		}
	})
	return pathScores
}

/* score branch u -> w (for w in subtree under sub); returns score, best w */
func (dp *DP) scoreU(u, sub, v *tree.Node, pathScores map[int]uint) (uint, int) {
	var bestScore uint
	var bestW int
	SubtreePreOrder(sub, func(w *tree.Node) {
		score := dp.scoreEdge(u, w, v, sub) + pathScores[w.Id()]
		if score >= bestScore {
			bestScore = score
			bestW = w.Id()
		}
	})
	return pathScores[u.Id()] + bestScore, bestW
}

func (dp *DP) scoreEdge(u, w, v, wSub *tree.Node) uint {
	score := uint(0)
	for _, q := range dp.TreeData.QuartetSet[v.Id()] {
		if dp.quartetScore(q, u, w, v, wSub) {
			score += 1
		}
	}
	return score
}

func (dp *DP) quartetScore(q *prep.Quartet, u, w, v, wSub *tree.Node) bool {
	bottom := -1
	bi := -1 // bottom index (0 - 3)
	for i, t := range q.Taxa {
		if dp.TreeData.Leafsets[w.Id()][t] && bottom == -1 {
			bottom = t
			bi = i
		} else if dp.TreeData.Leafsets[w.Id()][t] {
			return false
		}
	}
	lcaSet := make(map[uint]bool)
	for _, t := range q.Taxa {
		lcaSet[dp.TreeData.LCA[bottom][t]] = true
	}
	if len(lcaSet) != 4 {
		return false
	}
	neighbor := neighborTaxaQ(q, bi)
	var lcaDepths [4]int
	i := 0
	for k, v := range lcaSet {
		if v {
			d, err := dp.TreeData.IdToNodes[k].Depth()
			if err != nil {
				panic(err)
			}
			lcaDepths[i] = d
			i++
		}
	}
	nLeaves := len(dp.TreeData.Leafsets[0])
	maxW, minU, bestTaxa := -1, nLeaves, -1
	taxaInU := false
	for i, t := range q.Taxa {
		d, err := dp.TreeData.IdToNodes[t].Depth()
		if err != nil {
			panic(err)
		}
		if !taxaInU && dp.TreeData.Leafsets[wSub.Id()][t] && d > maxW {
			maxW = d
			bestTaxa = i
		} else if !dp.TreeData.Leafsets[wSub.Id()][t] && d < minU {
			taxaInU = true
			minU = d
			bestTaxa = i
		}
	}
	return (!taxaInU && q.Taxa[bestTaxa] == neighbor) || (taxaInU && q.Taxa[bestTaxa] == neighbor)
}

func min(vals [4]int) int {
	m := vals[0]
	for i := 1; i < 4; i++ {
		if vals[i] > m {
			m = vals[i]
		}
	}
	return m
}

/* return neighbor of taxa at index i in quartet */
func neighborTaxaQ(q *prep.Quartet, i int) int {
	b := (q.Topology >> (i - 1)) % 2
	for j := 0; j < 4; j++ {
		if (q.Topology>>(j-1))%2 == b {
			return j
		}
	}
	panic("invalid quartet!!")
}

func (dp *DP) traceback() [][2]int {
	return nil
}
