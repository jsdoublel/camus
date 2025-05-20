package infer

import "github.com/evolbioinfo/gotree/tree"

func SubtreePostOrder(cur *tree.Node, f func(cur, otherSubtree *tree.Node)) {
	if !cur.Tip() {
		children := make([]*tree.Node, 0)
		for _, n := range cur.Neigh() {
			if p, err := cur.Parent(); err != nil && err.Error() != "The node has no parent : May be the root?" {
				panic(err)
			} else if n != p {
				children = append(children, n)
			}
		}
		if len(children) != 2 {
			panic("tree is not binary")
		}
		subtreePostOrderHelper(children[0], children[1], f)
		subtreePostOrderHelper(children[1], children[0], f)
	}
}

func subtreePostOrderHelper(cur, otherSubtree *tree.Node, f func(cur, otherSubtree *tree.Node)) {
	for _, n := range cur.Neigh() {
		if p, err := cur.Parent(); err != nil {
			panic(err)
		} else if n != p {
			subtreePostOrderHelper(n, otherSubtree, f)
		}
	}
	f(cur, otherSubtree)
}

func SubtreePreOrder(cur *tree.Node, f func(cur *tree.Node)) {
	f(cur)
	for _, n := range cur.Neigh() {
		if p, err := cur.Parent(); err != nil && err.Error() != "The node has no parent : May be the root?" {
			panic(err)
		} else if n != p {
			SubtreePreOrder(n, f)
		}
	}
}

// returns best split between two lists, i.e., max l[i] + r[j] where i + j = k.
// returns err if k is too large.
func BestSplit(l, r []uint, k int) (int, int, error) {
	if len(l) == 0 || len(r) == 0 {
		panic("zero length lists not allowed")
	}
	if len(l)+len(r)-2 < k {
		return 0, 0, ErrNoValidSplit
	}
	bestScore := maxVal
	bestKL, bestKR := -1, -1
	for i := max(0, k-(len(r)-1)); i <= min(k, len(l)-1); i++ {
		if curScore := l[i] + r[k-i]; curScore > bestScore || bestScore == maxVal {
			bestScore = curScore
			bestKL, bestKR = i, k-i
		}
	}
	return bestKL, bestKR, nil
}

// Takes a value k, and four lists; returns a slice of indices
// idx = {idx1, idx2, ... } such that the sum of lists[0][idx1] + lists[1][idx2]
// + ... is maximized, and all indices add up to k.
// Returns an error if a valid split does not exist.
func FourWayBestSplit(lists [4][]uint, k int) (indices [4]int, err error) {
	combinedLen := 0
	for _, l := range lists {
		combinedLen += len(l)
		if len(l) == 0 {
			panic("zero length lists not allowed")
		}
	}
	if combinedLen-len(lists) < k {
		return [4]int{}, ErrNoValidSplit
	}
	sol1, scores1 := solveAllSplits(lists[0], lists[1], k)
	sol2, scores2 := solveAllSplits(lists[2], lists[3], k)
	idx1, idx2, err := BestSplit(scores1, scores2, k)
	if err != nil {
		panic("did not find expected split")
	}
	copy(indices[:2], sol1[idx1][:])
	copy(indices[2:], sol2[idx2][:])
	return
}

// helper for fourWayBestSplit
func solveAllSplits(list1, list2 []uint, k int) (solutions [][2]int, scores []uint) {
	solutions = make([][2]int, 0, k)
	scores = make([]uint, 0, k)
	for i := range k + 1 {
		l, r, err := BestSplit(list1, list2, i)
		if err != nil {
			break
		}
		solutions = append(solutions, [2]int{l, r})
		scores = append(scores, list1[l]+list2[r])
	}
	return
}
