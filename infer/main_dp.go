// Package implementing the CAMUS dynamic programming algorithm. We use the
// following naming convention for some variables throughout. v always
// represents the vertex corresponding to the current subproblem. At vertex v
// we consider adding a directed edge from u -> w.
package infer

import (
	"fmt"
	"log"

	"github.com/evolbioinfo/gotree/tree"

	gr "github.com/jsdoublel/camus/graphs"
	"github.com/jsdoublel/camus/prep"
)

type DP struct {
	DP           [][]uint      // score for each dp subproblem
	Branches     [][]gr.Branch // branch for each dp subproblem
	TreeData     *gr.TreeData  // preprocessed data for our constraint tree
	brScoreCache [][]uint      // caches branch scores so they are not recomputed (INDEXED [u][w] ONLY!!)
	NumNodes     int           // number of nodes
}

var noCache = ^uint(0)

// Runs CAMUS algorithm -- returns preprocessed tree data struct, quartet count stats, list of branches.
// Errors returned come from preprocessing (invalid inputs, etc.).
func CAMUS(tre *tree.Tree, geneTrees []*tree.Tree) (*gr.TreeData, [][]gr.Branch, error) {
	log.Print("beginning data preprocessing\n")
	td, err := prep.Preprocess(tre, geneTrees)
	if err != nil {
		return nil, nil, fmt.Errorf("preprocess error: %w", err)
	}
	log.Print("preprocessing finished, beginning dp algorithm\n")
	n := len(td.Tree.Nodes())
	dp := &DP{
		DP:           make([][]uint, n, n),
		Branches:     make([][]gr.Branch, n),
		TreeData:     td,
		brScoreCache: makeBrScoreCache(n),
		NumNodes:     n,
	}
	return td, dp.RunDP(), nil
}

func (dp *DP) RunDP() [][]gr.Branch {
	n, err := dp.TreeData.Tree.NbTips()
	if err != nil {
		panic(err)
	}
	count := 0
	totalDP := len(dp.DP) - n
	dp.TreeData.Tree.PostOrder(func(cur, prev *tree.Node, e *tree.Edge) (keep bool) {
		if !cur.Tip() { // TODO: fix no longer true: default value is 0, so we don't need to code a base case
			count++
			prep.LogEveryNPercent(count, 2, totalDP, fmt.Sprintf("processing dp cell %d of %d\n", count, totalDP))
			scores, branches := dp.score(cur)
			dp.DP[cur.Id()] = scores
			dp.Branches[cur.Id()] = branches
		} else {
			dp.DP[cur.Id()] = make([]uint, 1)
			dp.Branches[cur.Id()] = make([]gr.Branch, 1)
		}
		return true
	})
	numOptimal := len(dp.DP[dp.TreeData.Root.Id()]) - 1
	log.Printf("%d edges identified\n", numOptimal)
	log.Println("beginning traceback")
	result := make([][]gr.Branch, numOptimal, numOptimal)
	for m := range numOptimal + 1 {
		if m != 0 {
			finalScore := dp.DP[dp.TreeData.Tree.Root().Id()][m]
			log.Printf("dp scored %d at root with %d edges\n", finalScore, m)
			log.Printf("%f percent of quartets satisfied", 100*(float64(finalScore)/float64(dp.TreeData.TotalNumQuartets())))
			result[m-1] = dp.traceback(m)
		}
	}
	log.Println("done.")
	return result
}

func (dp *DP) score(v *tree.Node) ([]uint, []gr.Branch) {
	lID, rID := dp.TreeData.Children[v.Id()][0].Id(), dp.TreeData.Children[v.Id()][1].Id()
	scores := make([]uint, 1, dp.NumNodes) // choice of capacity is a bit arbitrary
	branches := make([]gr.Branch, 1, dp.NumNodes)
	scores[0] = dp.DP[lID][0] + dp.DP[rID][0]
	for m := 1; m < dp.NumNodes*dp.NumNodes; m++ { // N*N for safty: greater than possible number of branches
		score, branch := dp.scoreV(v, m)
		lM, rM := min(m, len(dp.DP[lID])-1), min(m, len(dp.DP[rID])-1)
		noEdgeScore := dp.DP[lID][lM] + dp.DP[rID][rM]
		if scores[m-1] >= score && scores[m-1] >= noEdgeScore { // converge
			break
		}
		if score > noEdgeScore {
			scores = append(scores, score)
			branches = append(branches, branch)
		} else {
			scores = append(scores, noEdgeScore)
			branches = append(branches, gr.Branch{})
		}
		if m == dp.NumNodes*dp.NumNodes-1 {
			panic("runaway loop")
		}
		if scores[m] <= scores[m-1] {
			panic("score did not strictly improve")
		}
		if len(scores) != len(branches) || len(scores) != m && len(scores) != m+1 {
			panic(fmt.Sprintf("scores list in weird state: m %d, len(scores) %d, len(branches) %d", m, len(scores), len(branches)))
		}
	}
	return scores, branches
}

// Calculates score for given top node v; returns score and best edge.
// m indicates that the edge being added is the m^th edge.
func (dp *DP) scoreV(v *tree.Node, m int) (uint, gr.Branch) {
	if m == 0 {
		panic("scoreV should never be called with zero m value")
	}
	prevM := m - 1
	bestScore := uint(0)
	var bestEdge gr.Branch
	bestCycleLen := 0
	vPathDPScores := dp.accumlateDPScores(v, prevM)
	if v != dp.TreeData.Root { // TODO: consider refactoring iterator to see if this case can be eliminated
		var wID int
		bestScore, wID = dp.scoreU(v, v, v, vPathDPScores, prevM)
		bestCycleLen = dp.cycleLen(v.Id(), wID)
		bestEdge = gr.Branch{IDs: [2]int{v.Id(), wID}}
	}
	SubtreePostOrder(v, func(cur, otherSubtree *tree.Node) { // might be slightly wrong in tie breaking: (larger cycles preferred in some cases)
		score, wID := dp.scoreU(cur, otherSubtree, v, vPathDPScores, prevM)
		if score > bestScore || (score == bestScore && dp.cycleLen(cur.Id(), wID) <= bestCycleLen) {
			bestScore, bestEdge = score, gr.Branch{IDs: [2]int{cur.Id(), wID}}
		}
	})
	return bestScore, bestEdge
}

// Add up dp scores along path from v ~> w
func (dp *DP) accumlateDPScores(v *tree.Node, prevM int) map[int]uint {
	pathScores := make(map[int]uint)
	SubtreePreOrder(v, func(cur *tree.Node) {
		if cur != v {
			if p, err := cur.Parent(); v != dp.TreeData.Root && err != nil {
				panic(err)
			} else if p != v && cur != dp.TreeData.Root {
				sibId := dp.TreeData.Sibling(cur).Id()
				mLookup := min(len(dp.DP[sibId])-1, prevM)
				pathScores[cur.Id()] = pathScores[p.Id()] + dp.DP[sibId][mLookup]
			}
		}
	})
	return pathScores
}

func (dp *DP) cycleLen(u, w int) int {
	v := dp.TreeData.LCA(u, w)
	length := (dp.TreeData.Depths[u] - dp.TreeData.Depths[v]) + (dp.TreeData.Depths[w] - dp.TreeData.Depths[v]) + 1
	if v == u { // we have to account for the edge above v that our new edge is anchored to
		length += 1
	}
	return length
}

// Score branch u -> w (for w in subtree under sub); returns score, best w
func (dp *DP) scoreU(u, sub, v *tree.Node, pathScores map[int]uint, prevM int) (uint, int) {
	var bestScore uint
	var bestW int
	SubtreePreOrder(sub, func(w *tree.Node) {
		if u != w {
			edgeScore := dp.brScoreCache[u.Id()][w.Id()]
			if edgeScore == noCache {
				edgeScore = dp.scoreEdge(u, w, v, sub)
				dp.brScoreCache[u.Id()][w.Id()] = edgeScore
			}
			mLookup := min(len(dp.DP[w.Id()])-1, prevM)
			score := edgeScore + pathScores[w.Id()] + dp.DP[w.Id()][mLookup]
			if score > bestScore {
				bestScore = score
				bestW = w.Id()
			}
		}
	})
	score := pathScores[u.Id()] + bestScore
	if u != v {
		mLookup := min(len(dp.DP[u.Id()])-1, prevM)
		score += dp.DP[u.Id()][mLookup]
	}
	return score, bestW
}

func (dp *DP) scoreEdge(u, w, v, wSub *tree.Node) uint {
	score := uint(0)
	for _, q := range dp.TreeData.Quartets(v.Id()) {
		if QuartetScore(q, u, w, v, wSub, dp.TreeData) == gr.Qeq {
			score += dp.TreeData.NumQuartet(q)
		}
	}
	return score
}

func QuartetScore(q *gr.Quartet, u, w, v, wSub *tree.Node, td *gr.TreeData) int {
	bottom, bi, unique := uniqueTaxaBelowNodeFromQ(w, q, td)
	if !unique || bottom == -1 {
		return gr.Qdiff
	}
	cycleNodes := make(map[int]bool)
	taxaToLCA := make(map[int]int) // tip index -> lca
	for _, t := range q.Taxa {
		tID := td.NodeID(t)
		var lca int
		if !td.InLeafset(v.Id(), t) {
			lca = 0
		} else if td.InLeafset(wSub.Id(), t) || td.InLeafset(u.Id(), bottom) {
			lca = td.LCA(w.Id(), tID)
		} else {
			lca = td.LCA(u.Id(), tID)
		}
		cycleNodes[lca] = true
		taxaToLCA[t] = lca
	}
	if len(cycleNodes) != 4 {
		return gr.Qdiff
	}
	neighbor := neighborTaxaQ(q, bi)
	lcaDepths := make(map[int]int) // node ID -> depth
	for k, v := range cycleNodes {
		if v {
			lcaDepths[k] = td.Depths[k]
		}
	}
	nLeaves := td.NLeaves
	minW, maxU, bestTaxa := nLeaves, -1, -1
	taxaInU := false
	for _, t := range q.Taxa {
		d := lcaDepths[taxaToLCA[t]]
		if !taxaInU && (td.InLeafset(wSub.Id(), t) && d < minW) {
			minW = d
			bestTaxa = t
		} else if !td.InLeafset(wSub.Id(), t) && d > maxU {
			taxaInU = true
			maxU = d
			bestTaxa = t
		}
	}
	if bestTaxa == neighbor {
		return gr.Qeq
	} else {
		return gr.Qneq
	}
}

// Returns -1 for both id and index if no taxa is found, true if taxa is unique (or there isn't a taxa)
func uniqueTaxaBelowNodeFromQ(n *tree.Node, q *gr.Quartet, td *gr.TreeData) (int, int, bool) {
	taxaID, taxaIndex := -1, -1
	for i, t := range q.Taxa {
		if td.InLeafset(n.Id(), t) && taxaID == -1 {
			taxaID, taxaIndex = t, i
		} else if td.InLeafset(n.Id(), t) {
			return taxaID, taxaIndex, false
		}
	}
	return taxaID, taxaIndex, true
}

// Return neighbor of taxa at index i in quartet
func neighborTaxaQ(q *gr.Quartet, i int) int {
	b := (q.Topology >> i) % 2
	for j := range 4 {
		if j != i && (q.Topology>>j)%2 == b {
			return q.Taxa[j]
		}
	}
	panic("invalid quartet or bad i")
}

func (dp *DP) traceback(m int) []gr.Branch {
	return dp.tracebackRecursive(dp.TreeData.Root, m)
}

func (dp *DP) tracebackRecursive(curNode *tree.Node, m int) []gr.Branch {
	if !dp.TreeData.IdToNodes[curNode.Id()].Tip() {
		mLookup := min(len(dp.Branches[curNode.Id()])-1, m)
		curBranch := dp.Branches[curNode.Id()][mLookup]
		if curBranch.Empty() {
			return append(dp.tracebackRecursive(dp.TreeData.Children[curNode.Id()][0], m),
				dp.tracebackRecursive(dp.TreeData.Children[curNode.Id()][1], m)...)
		} else {
			prevM := m - 1
			u, w := dp.TreeData.IdToNodes[curBranch.IDs[0]], dp.TreeData.IdToNodes[curBranch.IDs[1]]
			traceback := []gr.Branch{curBranch}
			if u != curNode {
				traceback = append(traceback, dp.tracebackRecursive(u, prevM)...)
				traceback = append(traceback, dp.tracePath(u, curNode, prevM)...)
			}
			if w == curNode {
				panic("w should not be current node")
			}
			traceback = append(traceback, dp.tracebackRecursive(w, prevM)...)
			traceback = append(traceback, dp.tracePath(w, curNode, prevM)...)
			return traceback
		}
	}
	return []gr.Branch{}
}

func (dp *DP) tracePath(start, end *tree.Node, prevM int) []gr.Branch {
	if start == end {
		panic("in tracePath start should not equal end!")
	}
	cur := start
	p, err := start.Parent()
	if err != nil {
		panic(err)
	}
	trace := make([]gr.Branch, 0)
	for p != end {
		trace = append(trace, dp.tracebackRecursive(dp.TreeData.Sibling(cur), prevM)...)
		var err error
		cur = p
		p, err = cur.Parent()
		if err != nil {
			panic(err)
		}
	}
	return trace
}

func makeBrScoreCache(n int) [][]uint {
	result := make([][]uint, n, n)
	for i := range result {
		result[i] = make([]uint, n, n)
		for j := range result[i] {
			result[i][j] = noCache
		}
	}
	return result
}
