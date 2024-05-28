package alg

import (
	"fmt"
	"os"

	"github.com/evolbioinfo/gotree/tree"

	"camus/prep"
)

func CAMUS(tre *tree.Tree, rawQuartets []*tree.Tree) ([][2]int, error) {
	fmt.Fprint(os.Stderr, "beginning data preprocessing")
	quartets, td, err := prep.Preprocess(tre, rawQuartets)
	if err != nil {
		return nil, fmt.Errorf("Preprocess error: %w", err)
	}
	fmt.Fprint(os.Stderr, "preprocessing finished, beginning CAMUS")
	n := len(td.Tree.Nodes())
	dp := make([]uint, n, n)      // score for each dp subproblem
	branches := make([][2]int, n) // branch for each dp subproblem
	backtrace := make([][]int, n) // prev subproblems for each dp subproblem
	td.Tree.PostOrder(func(cur, prev *tree.Node, e *tree.Edge) (keep bool) {
		if !cur.Tip() { // default value is 0, so we don't need to code a base case
			lID, rID := td.Children[cur.Id()][0].Id(), td.Children[cur.Id()][1].Id()
			score, branch, back := score(cur, dp, td, quartets)
			noEdgeScore := dp[lID] + dp[rID]
			if score > noEdgeScore {
				dp[cur.Id()] = score
				branches[cur.Id()] = branch
				backtrace[cur.Id()] = back
			} else {
				dp[cur.Id()] = noEdgeScore
				backtrace[cur.Id()] = []int{lID, rID}
			}
		}
		return true
	})
	result := traceback(branches, backtrace)
	return result, nil
}

func score(alpha *tree.Node, dp []uint, td *prep.TreeData, qSet []*prep.Quartet) (uint, [2]int, []int) {
	return 0, [2]int{0, 0}, []int{}
}

func traceback(branches [][2]int, backtrace [][]int) [][2]int {
	return nil
}
