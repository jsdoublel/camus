package alg

import (
	"fmt"
	"log"

	"github.com/evolbioinfo/gotree/tree"

	"camus/prep"
)

type DP struct {
	DP       []uint         // score for each dp subproblem
	Branches [][2]int       // branch for each dp subproblem
	TreeData *prep.TreeData // preprocessed data for our constraint tree
}

func CAMUS(tre *tree.Tree, geneTrees []*tree.Tree) (*prep.TreeData, [][2]int, error) {
	log.Print("beginning data preprocessing\n")
	td, err := prep.Preprocess(tre, geneTrees)
	if err != nil {
		return nil, nil, fmt.Errorf("preprocess error: %w", err)
	}
	log.Print("preprocessing finished, beginning CAMUS\n")
	n := len(td.Tree.Nodes())
	dp := &DP{DP: make([]uint, n, n), Branches: make([][2]int, n), TreeData: td}
	return td, dp.RunDP(), nil
}

func (dp *DP) RunDP() [][2]int {
	count := 0
	dp.TreeData.Tree.PostOrder(func(cur, prev *tree.Node, e *tree.Edge) (keep bool) {
		if !cur.Tip() { // default value is 0, so we don't need to code a base case
			count++
			log.Printf("processing dp cell %d of %d\n", count, len(dp.DP)-len(dp.TreeData.Tree.Tips()))
			lID, rID := dp.TreeData.Children[cur.Id()][0].Id(), dp.TreeData.Children[cur.Id()][1].Id()
			score, branch := dp.score(cur)
			noEdgeScore := dp.DP[lID] + dp.DP[rID]
			if score > noEdgeScore {
				dp.DP[cur.Id()] = score
				dp.Branches[cur.Id()] = branch
			} else {
				dp.DP[cur.Id()] = noEdgeScore
			}
		}
		return true
	})
	log.Printf("dp scored %d at root\n", dp.DP[dp.TreeData.Tree.Root().Id()])
	result := dp.traceback()
	log.Printf("%d edges identified\n", len(result))
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
		if cur != v {
			if p, err := cur.Parent(); v != dp.TreeData.Root && err != nil {
				panic(err)
			} else if p != v && cur != dp.TreeData.Root {
				pathScores[cur.Id()] = pathScores[p.Id()] + dp.DP[dp.TreeData.Sibling(cur).Id()]
			}
		}
	})
	return pathScores
}

/* score branch u -> w (for w in subtree under sub); returns score, best w */
func (dp *DP) scoreU(u, sub, v *tree.Node, pathScores map[int]uint) (uint, int) {
	var bestScore uint
	var bestW int
	SubtreePreOrder(sub, func(w *tree.Node) {
		if u != w {
			score := dp.scoreEdge(u, w, v, sub) + pathScores[w.Id()] + dp.DP[w.Id()]
			if score >= bestScore {
				bestScore = score
				bestW = w.Id()
			}
		}
	})
	return pathScores[u.Id()] + bestScore, bestW
}

func (dp *DP) scoreEdge(u, w, v, wSub *tree.Node) uint {
	score := uint(0)
	for _, q := range dp.TreeData.QuartetSet[v.Id()] {
		if dp.quartetScore(q, u, w, v, wSub) {
			score += dp.TreeData.QuartetCounts[*q]
		}
	}
	return score
}

func (dp *DP) quartetScore(q *prep.Quartet, u, w, v, wSub *tree.Node) bool {
	bottom, bi, unique := dp.uniqueTaxaBelowNodeFromQ(w, q)
	if !unique || bottom == -1 {
		return false
	}
	cycleNodes := make(map[int]bool)
	taxaToLCA := make(map[int]int) // tip index -> lca
	for _, t := range q.Taxa {
		tID := dp.TreeData.NodeID(t)
		var lca int
		if !dp.TreeData.InLeafset(v.Id(), t) {
			lca = 0
		} else if dp.TreeData.InLeafset(wSub.Id(), t) || dp.TreeData.InLeafset(u.Id(), bottom) {
			lca = dp.TreeData.LCA(w.Id(), tID)
		} else {
			lca = dp.TreeData.LCA(u.Id(), tID)
		}
		cycleNodes[lca] = true
		taxaToLCA[t] = lca
	}
	if len(cycleNodes) != 4 {
		return false
	}
	neighbor := neighborTaxaQ(q, bi)
	lcaDepths := make(map[int]int) // node ID -> depth
	for k, v := range cycleNodes {
		if v {
			lcaDepths[k] = dp.TreeData.Depths[k]
		}
	}
	nLeaves := dp.TreeData.NLeaves()
	minW, maxU, bestTaxa := nLeaves, -1, -1
	taxaInU := false
	for _, t := range q.Taxa {
		d := lcaDepths[taxaToLCA[t]]
		if !taxaInU && (dp.TreeData.InLeafset(wSub.Id(), t) && d < minW) {
			minW = d
			bestTaxa = t
		} else if !dp.TreeData.InLeafset(wSub.Id(), t) && d > maxU {
			taxaInU = true
			maxU = d
			bestTaxa = t
		}
	}
	return bestTaxa == neighbor
}

/* returns -1 for both id and index if no taxa is found, true if taxa is unique (or there isn't a taxa) */
func (dp *DP) uniqueTaxaBelowNodeFromQ(n *tree.Node, q *prep.Quartet) (int, int, bool) {
	taxaID, taxaIndex := -1, -1
	for i, t := range q.Taxa {
		if dp.TreeData.InLeafset(n.Id(), t) && taxaID == -1 {
			taxaID, taxaIndex = t, i
		} else if dp.TreeData.InLeafset(n.Id(), t) {
			return taxaID, taxaIndex, false
		}
	}
	return taxaID, taxaIndex, true
}

/* return neighbor of taxa at index i in quartet */
func neighborTaxaQ(q *prep.Quartet, i int) int {
	b := (q.Topology >> i) % 2
	for j := range 4 {
		if j != i && (q.Topology>>j)%2 == b {
			return q.Taxa[j]
		}
	}
	panic("invalid quartet or bad i")
}

func (dp *DP) traceback() [][2]int {
	return dp.tracebackRecursive(dp.TreeData.Root)
}

func (dp *DP) tracebackRecursive(curNode *tree.Node) [][2]int {
	if !dp.TreeData.IdToNodes[curNode.Id()].Tip() {
		curBranch := dp.Branches[curNode.Id()]
		if curBranch == [2]int{0, 0} { // there is no edge
			return append(dp.tracebackRecursive(dp.TreeData.Children[curNode.Id()][0]),
				dp.tracebackRecursive(dp.TreeData.Children[curNode.Id()][1])...)
		} else {
			u, w := dp.TreeData.IdToNodes[curBranch[0]], dp.TreeData.IdToNodes[curBranch[1]]
			traceback := [][2]int{curBranch}
			if u != curNode {
				traceback = append(traceback, dp.tracebackRecursive(u)...)
				traceback = append(traceback, dp.tracePath(u, curNode)...)
			}
			if w == curNode {
				panic("w should not be current node")
			}
			traceback = append(traceback, dp.tracebackRecursive(w)...)
			traceback = append(traceback, dp.tracePath(w, curNode)...)
			return traceback
		}
	}
	return [][2]int{}
}

func (dp *DP) tracePath(start, end *tree.Node) [][2]int {
	if start == end {
		panic("in tracePath start should not equal end!")
	}
	cur := start
	p, err := start.Parent()
	if err != nil {
		panic(err)
	}
	trace := make([][2]int, 0)
	for p != end {
		trace = append(trace, dp.tracebackRecursive(dp.TreeData.Sibling(cur))...)
		var err error
		cur = p
		p, err = cur.Parent()
		if err != nil {
			panic(err)
		}
	}
	return trace
}
