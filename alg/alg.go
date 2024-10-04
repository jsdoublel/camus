package alg

import (
	"fmt"
	"os"

	"github.com/evolbioinfo/gotree/tree"

	"camus/prep"
)

func CAMUS(tre *tree.Tree, geneTrees []*tree.Tree) ([][2]int, error) {
	fmt.Fprint(os.Stderr, "beginning data preprocessing")
	td, err := prep.Preprocess(tre, geneTrees)
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
			score, branch, back := score(cur, dp, td)
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

func score(alpha *tree.Node, dp []uint, td *prep.TreeData) (uint, [2]int, []int) {
	// loop over all possible values of u for a given alpha
	//		scoreU(u, ...) + dp lookup on path from v to u
	return 0, [2]int{0, 0}, []int{}
}

/* score branch u -> w (for w in subtree under sub); returns score, best w */
func scoreU(u, sub *tree.Node, qSet []*prep.Quartet, td *prep.TreeData) (uint, int) {
	// solve topology for every quartet in the quartet set
	//		loop over all leaves in subtree
	//		score edge directly above leaf
	return 0, 0
}

func traceback(branches [][2]int, backtrace [][]int) [][2]int {
	return nil
}
