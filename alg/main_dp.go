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
	return dp.RunDP()
}

func (dp *DP) RunDP() ([][2]int, error) {
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
	return result, nil
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
	SubtreePreOrder(sub, func(cur *tree.Node) {
		score := dp.scoreEdge(u, cur, v) + pathScores[cur.Id()]
		if score >= bestScore {
			bestScore = score
			bestW = cur.Id()
		}
	})
	return pathScores[u.Id()] + bestScore, bestW
}

func (dp *DP) scoreEdge(u, w, v *tree.Node) uint {
	score := uint(0)
	for _, q := range dp.TreeData.QuartetSet[v.Id()] {
		if dp.quartetScore(q, u, w, v) {
			score += 1
		}
	}
	return score
}

func (dp *DP) quartetScore(q *prep.Quartet, u, w, v *tree.Node) bool {
	numLeaves := uint(len(dp.TreeData.Leafsets[0]))
	bottom := numLeaves
	bIdx := -1
	for i, t := range q.Taxa {
		if dp.TreeData.Leafsets[w.Id()][t] && bottom == numLeaves {
			bottom = t
			bIdx = i
		} else if dp.TreeData.Leafsets[w.Id()][t] {
			return false
		}
	}
	// TODO: check that all taxa are on separate subtrees coming off of the cycle
	neighbor := neighborTaxaQ(q, bIdx)
	vl, vr := dp.TreeData.Children[v.Id()][0], dp.TreeData.Children[v.Id()][1]
	return false
}

func (dp *DP) traceback() [][2]int {
	return nil
}
