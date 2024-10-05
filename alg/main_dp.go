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
			score, branch, back := dp.score(cur)
			noEdgeScore := dp.DP[lID] + dp.DP[rID]
			if score > noEdgeScore {
				dp.DP[cur.Id()] = score
				dp.Branches[cur.Id()] = branch
				dp.Backtrace[cur.Id()] = back
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
func (dp *DP) score(v *tree.Node) (uint, [2]int, []int) {
	var bestScore uint
	var bestEdge [2]int
	var bestTrace []int
	vPathDPScores := dp.accumlateDPScores(v) // TODO: think if there's a better way to update this
	if v != dp.TreeData.Root {
		var wID int
		bestScore, wID, bestTrace = dp.scoreU(v, v, v, vPathDPScores)
		bestEdge = [2]int{v.Id(), wID}
	}
	SubtreePostOrder(v, func(cur, otherSubtree *tree.Node) {
		score, wID, back := dp.scoreU(cur, otherSubtree, v, vPathDPScores)
		if score >= bestScore {
			bestScore, bestEdge, bestTrace = score, [2]int{cur.Id(), wID}, back
		}
	})
	return bestScore, bestEdge, bestTrace
}

// TODO: add backtrace to this part
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
func (dp *DP) scoreU(u, sub, v *tree.Node, pathScores map[int]uint) (uint, int, []int) {
	// qSet := dp.TreeData.QuartetSets[v.Id()]
	// // iterate through possible w
	// var bestScore uint
	// var bestW int
	// SubtreePreOrder(sub, func(cur *tree.Node) {
	// })
	return 0, 0, []int{}
}

func (dp *DP) scoreEdge(u, w *tree.Node) uint {
	return 0
}

func (dp *DP) traceback() [][2]int {
	return nil
}
