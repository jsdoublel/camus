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
	TreeData  *prep.TreeData // preprocessed data for our tree
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
	td.Tree.PostOrder(func(cur, prev *tree.Node, e *tree.Edge) (keep bool) {
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
	// loop over all possible values of u for a given v
	//		scoreU(u, ...) + dp lookup on path from v to u
	scores := make([]uint, 0)
	edges := make([][2]int, 0)
	backtraces := make([][]int, 0)
	if v != dp.TreeData.Root {
		score, w, back := dp.scoreU(v, v, v)
		scores = append(scores, score)
		edges = append(edges, [2]int{v.Id(), w})
		backtraces = append(backtraces, back)
	}
	PostOrderU(v, func(cur, otherSubtree *tree.Node) {
		// do dp lookups along path between v -> u
		score, w, back := dp.scoreU(cur, otherSubtree, v)
		scores = append(scores, score)
		edges = append(edges, [2]int{v.Id(), w})
		backtraces = append(backtraces, back)
	})
	maxScore, maxIdx := 0, 0
	for i, s := range scores {
		if s >= uint(maxScore) {
			maxScore = int(s)
			maxIdx = i
		}
	}
	return uint(maxScore), edges[maxIdx], backtraces[maxIdx]
}

/* score branch u -> w (for w in subtree under sub); returns score, best w */
func (dp *DP) scoreU(u, sub, v *tree.Node) (uint, int, []int) {
	qSet := dp.TreeData.QuartetSets[v.Id()]
	// iterate through possible w
	scores := make([]int, 0)
	backtraces := make([][]int, 0)
	SubtreePreOrder(sub, func(cur *tree.Node) {

	})
	return 0, 0, []int{}
}

func PostOrderU(cur *tree.Node, f func(cur, otherSubtree *tree.Node, lookupTotal uint)) {
	if !cur.Tip() {
		children := make([]*tree.Node, 0)
		for _, n := range cur.Neigh() {
			if p, err := cur.Parent(); err != nil {
				panic(err)
			} else if n != p {
				children = append(children, n)
			}
		}
		if len(children) != 2 {
			panic("tree is not binary")
		}
		subtreePostOrderHelper(children[0], children[1], 0, f)
		subtreePostOrderHelper(children[1], children[0], 0, f)
	} else {
	}

}

func subtreePostOrderHelper(cur, otherSubtree *tree.Node, lookupTotal uint, f func(cur, otherSubtree *tree.Node, lookupTotal uint)) {
	for _, n := range cur.Neigh() {
		if p, err := cur.Parent(); err != nil {
			panic(err)
		} else if n != p {
			subtreePostOrderHelper(n, otherSubtree, lookupTotal, f)
		}
	}
	f(cur, otherSubtree)
}

func SubtreePreOrder(cur *tree.Node, f func(cur *tree.Node)) {
	f(cur)
	for _, n := range cur.Neigh() {
		if p, err := cur.Parent(); err != nil {
			panic(err)
		} else if n != p {
			SubtreePreOrder(n, f)
		}
	}
}

func scoreW() {
}

// func solveTopology(q *prep.Quartet, u, taxa *tree.Node, td *prep.TreeData) bool {
// 	return false
// }

func (dp *DP) traceback() [][2]int {
	return nil
}
