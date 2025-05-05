// Package implementing the CAMUS dynamic programming algorithm. We use the
// following naming convention for some variables throughout. v always
// represents the vertex corresponding to the current subproblem. At vertex v
// we consider adding a directed edge from u -> w.
// k refers to the maximum number of edges we are considering at a given subproblem.
package infer

import (
	"errors"
	"fmt"
	"log"

	"github.com/evolbioinfo/gotree/tree"

	gr "github.com/jsdoublel/camus/graphs"
	"github.com/jsdoublel/camus/prep"
)

const maxVal = ^uint(0)

var ErrNoValidSplit = errors.New("no valid split")

type DP struct {
	DP           [][]uint     // score for each dp subproblem
	Traceback    [][]trace    // traceback for each dp subproblem
	TreeData     *gr.TreeData // preprocessed data for our constraint tree
	brScoreCache [][]uint     // caches branch scores so they are not recomputed (INDEXED [u][w] ONLY!!)
	NumNodes     int          // number of nodes
}

// stores accumlated dp scores for v ~> w paths
type pathDPScores struct {
	scores     []uint
	traceNodes []*cycleTraceNode
}

func (ps *pathDPScores) set(i int, score uint, traceNode cycleTraceNode) {
	ps.scores[i] = score
	ps.traceNodes[i] = &traceNode
}

func (ps pathDPScores) get(i int) (uint, *cycleTraceNode) {
	return ps.scores[i], ps.traceNodes[i]
}

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
		DP:           make([][]uint, n),
		Traceback:    make([][]trace, n),
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
	dp.TreeData.Tree.PostOrder(func(v, prev *tree.Node, e *tree.Edge) (keep bool) {
		if !v.Tip() {
			count++
			prep.LogEveryNPercent(count, 2, totalDP, fmt.Sprintf("processing dp cell %d of %d\n", count, totalDP))
			scores, edgeTrace := dp.solve(v)
			dp.DP[v.Id()] = scores
			dp.Traceback[v.Id()] = edgeTrace
		} else {
			dp.DP[v.Id()] = make([]uint, 1)
			dp.Traceback[v.Id()] = make([]trace, 1, dp.NumNodes)
			dp.Traceback[v.Id()][0] = noCycleTrace{}
		}
		return true
	})
	numOptimal := len(dp.DP[dp.TreeData.Root.Id()]) - 1
	log.Printf("%d edges identified\n", numOptimal)
	log.Println("beginning traceback")
	result := make([][]gr.Branch, numOptimal)
	for k := range numOptimal + 1 {
		if k != 0 {
			finalScore := dp.DP[dp.TreeData.Tree.Root().Id()][k]
			log.Printf("dp scored %d at root with %d edges\n", finalScore, k)
			log.Printf("%f percent of quartets satisfied", 100*(float64(finalScore)/float64(dp.TreeData.TotalNumQuartets())))
			result[k-1] = dp.traceback(k)
		}
	}
	log.Println("done.")
	return result
}

func (dp *DP) solve(v *tree.Node) ([]uint, []trace) {
	lID, rID := dp.TreeData.Children[v.Id()][0].Id(), dp.TreeData.Children[v.Id()][1].Id()
	scores := make([]uint, 1, dp.NumNodes) // choice of capacity is a bit arbitrary
	traces := make([]trace, 1, dp.NumNodes)
	scores[0] = dp.DP[lID][0] + dp.DP[rID][0]
	traces[0] = noCycleTrace{[2]*trace{&dp.Traceback[lID][0], &dp.Traceback[rID][0]}}
	for k := 1; ; k++ {
		score, edgeTrace := dp.scoreV(v, k)
		lK, rK, err := bestSplit(dp.DP[lID], dp.DP[rID], k)
		noEdgeScore := dp.DP[lID][lK] + dp.DP[rID][rK]
		if scores[k-1] >= score && (err != nil || scores[k-1] >= noEdgeScore) {
			break
		}
		if score > noEdgeScore {
			scores = append(scores, score)
			traces = append(traces, edgeTrace)
		} else if err == nil {
			scores = append(scores, noEdgeScore)
			traces = append(traces, &noCycleTrace{
				prevs: [2]*trace{&dp.Traceback[lID][lK], &dp.Traceback[rID][rK]},
			})
		}
		if k == dp.NumNodes*dp.NumNodes {
			panic("runaway loop")
		}
		if scores[k] <= scores[k-1] {
			panic("score did not strictly improve")
		}
		if len(scores) != len(traces) || len(scores) != k && len(scores) != k+1 {
			panic(fmt.Sprintf("scores list in weird state: k %d, len(scores) %d, len(branches) %d", k, len(scores), len(traces)))
		}
	}
	return scores, traces
}

// returns best split between two lists, i.e., max l[i] + r[j] where i + j = k.
// returns err if k is too large
func bestSplit(l, r []uint, k int) (int, int, error) {
	if len(l) == 0 || len(r) == 0 {
		panic("zero length lists not allowed")
	}
	if len(l)+len(r)-2 < k {
		return 0, 0, ErrNoValidSplit
	}
	bestScore := maxVal
	bestKL, bestKR := -1, -1
	for i := range min(k+1, len(l)) {
		if k-i >= len(r) {
			continue
		}
		if curScore := l[i] + r[k-i]; curScore > bestScore || bestScore == maxVal {
			bestScore = curScore
			bestKL, bestKR = i, k-i
		}
	}
	return bestKL, bestKR, nil
}

// Calculates score for given top node v; returns score and best edge.
// k indicates that the edge being added is the k^th edge.
func (dp *DP) scoreV(v *tree.Node, k int) (uint, trace) {
	if k == 0 {
		panic("scoreV should never be called with zero k value")
	}
	prevK := k - 1
	bestScore := uint(0)
	bestCycleTrace := cycleTrace{}
	bestCycleLen := 0
	vPathDPScores := dp.accumlateDPScores(v, prevK)
	if v != dp.TreeData.Root {
		var wID int
		bestScore, wID, bestCycleTrace = dp.scoreU(v, v, v, vPathDPScores, prevK)
		bestCycleLen = dp.cycleLen(v.Id(), wID)
		// curCycleTrace.branch, curCycleTrace.prevs = gr.Branch{IDs: [2]int{v.Id(), wID}}, prevs
	}
	SubtreePostOrder(v, func(u, otherSubtree *tree.Node) {
		score, wID, curCycleTrace := dp.scoreU(u, otherSubtree, v, vPathDPScores, prevK)
		if score > bestScore || (score == bestScore && dp.cycleLen(u.Id(), wID) <= bestCycleLen) {
			bestScore = score
			bestCycleTrace = curCycleTrace
			// curCycleTrace.branch, curCycleTrace.prevs = gr.Branch{IDs: [2]int{u.Id(), wID}}, prevs
		}
	})
	return bestScore, bestCycleTrace
}

// Add up dp scores along path from v ~> w
func (dp *DP) accumlateDPScores(v *tree.Node, prevK int) pathDPScores {
	pathScores := pathDPScores{
		make([]uint, dp.NumNodes),
		make([]*cycleTraceNode, dp.NumNodes),
	}
	SubtreePreOrder(v, func(cur *tree.Node) {
		if cur != v {
			if p, err := cur.Parent(); v != dp.TreeData.Root && err != nil {
				panic(err)
			} else if p != v && cur != dp.TreeData.Root {
				sibId := dp.TreeData.Sibling(cur).Id()
				kLookup := min(len(dp.DP[sibId])-1, prevK)
				pScore, pTrace := pathScores.get(p.Id())
				pathScores.set(
					cur.Id(),
					pScore+dp.DP[sibId][kLookup],
					cycleTraceNode{p: pTrace, sib: &dp.Traceback[sibId][kLookup]},
				)
				// pathScores[cur.Id()] = pathScores[p.Id()] + dp.DP[sibId][kLookup]
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
func (dp *DP) scoreU(u, sub, v *tree.Node, pathScores pathDPScores, prevK int) (uint, int, cycleTrace) {
	bestScore := maxVal
	var bestW int
	var bestWKLookup int
	var bestWTrace *cycleTraceNode
	SubtreePreOrder(sub, func(w *tree.Node) {
		if u != w { // TODO: is this check necessary?
			edgeScore := dp.brScoreCache[u.Id()][w.Id()]
			if edgeScore == maxVal {
				edgeScore = dp.scoreEdge(u, w, v, sub)
				dp.brScoreCache[u.Id()][w.Id()] = edgeScore
			}
			kLookup := min(len(dp.DP[w.Id()])-1, prevK)
			wScore, wTrace := pathScores.get(w.Id())
			score := edgeScore + wScore + dp.DP[w.Id()][kLookup]
			if score > bestScore || bestScore == maxVal {
				bestScore = score
				bestW = w.Id()
				bestWKLookup = kLookup
				bestWTrace = wTrace
			}
		}
	})
	uScore, uTrace := pathScores.get(u.Id())
	score := uScore + bestScore
	if u != v {
		mLookup := min(len(dp.DP[u.Id()])-1, prevK)
		score += dp.DP[u.Id()][mLookup]
	}
	return score, bestW, cycleTrace{
		pathW:  bestWTrace,
		pathU:  uTrace,
		wTrace: &dp.Traceback[bestW][bestWKLookup],
		branch: gr.Branch{IDs: [2]int{u.Id(), bestW}},
	}
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

func (dp DP) traceback(k int) []gr.Branch {
	return dp.Traceback[dp.TreeData.Root.Id()][k].traceback()
}

func makeBrScoreCache(n int) [][]uint {
	result := make([][]uint, n)
	for i := range result {
		result[i] = make([]uint, n)
		for j := range result[i] {
			result[i][j] = maxVal
		}
	}
	return result
}
