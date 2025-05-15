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

// Stores main dp algorithm data
type DP struct {
	DP           [][]uint     // score for each dp subproblem (DP[v][k])
	Traceback    [][]trace    // traceback for each dp subproblem (Traceback[v][k])
	TreeData     *gr.TreeData // preprocessed data for our constraint tree
	NumNodes     int          // number of nodes
	brScoreCache [][]uint     // caches branch scores so they are not recomputed (INDEXED [u][w] ONLY!!)
}

// Stores DP info for lookups corresponding to a given vertex v
type cycleDPLookups struct {
	v          *tree.Node
	scores     [][]uint            // score for each path (scores[w][k]); unique struct exists for each v
	traceNodes [][]*cycleTraceNode // backtrace for each path (traceNodes[w][k])
}

func (cl *cycleDPLookups) solveK(prevK int, dp *DP) {
	SubtreePreOrder(cl.v, func(cur *tree.Node) {
		if prevK == 0 {
			cl.scores[cur.Id()] = make([]uint, 0)
			cl.traceNodes[cur.Id()] = make([]*cycleTraceNode, 0)
		}
		if cur == cl.v { // don't want to look at parent of root/v
			return
		}
		p, err := cur.Parent()
		if cl.v != dp.TreeData.Root && err != nil {
			panic(err)
		} else if p == cl.v { // if parent is v, then sibling node of cur is also in the cycle
			return
		}
		sibId := dp.TreeData.Sibling(cur).Id()
		pScores, pTraces := cl.scores[p.Id()], cl.traceNodes[p.Id()]
		pK, sibK, err := BestSplit(pScores, dp.DP[sibId], prevK)
		if err != nil {
			return
		}
		cl.set(
			cur.Id(),
			prevK,
			pScores[pK]+dp.DP[sibId][sibK],
			cycleTraceNode{p: pTraces[pK], sib: &dp.Traceback[sibId][sibK]},
		)
	})
}

func (cl *cycleDPLookups) set(i, k int, score uint, traceNode cycleTraceNode) {
	cl.scores[i][k] = score // TODO: these probably need to be appends
	cl.traceNodes[i][k] = &traceNode
}

func (cl *cycleDPLookups) get(i, k int) (uint, *cycleTraceNode) {
	return cl.scores[i][k], cl.traceNodes[i][k]
}

// Runs CAMUS algorithm -- returns preprocessed tree data struct, quartet count stats, list of branches.
// Errors returned come from preprocessing (invalid inputs, etc.).
func CAMUS(tre *tree.Tree, geneTrees []*tree.Tree) (*gr.TreeData, [][]gr.Branch, error) {
	log.Println("beginning data preprocessing")
	td, err := prep.Preprocess(tre, geneTrees)
	if err != nil {
		return nil, nil, fmt.Errorf("preprocess error: %w", err)
	}
	log.Println("preprocessing finished, beginning dp algorithm")
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
	vCycleDP := cycleDPLookups{
		v:          v,
		scores:     make([][]uint, dp.NumNodes),
		traceNodes: make([][]*cycleTraceNode, dp.NumNodes),
	}
	for k := 1; ; k++ {
		score, edgeTrace, errAdd := dp.scoreV(v, k, &vCycleDP) // TODO: refactor into a function to make only a single error
		lK, rK, errNoAdd := BestSplit(dp.DP[lID], dp.DP[rID], k)
		noEdgeScore := dp.DP[lID][lK] + dp.DP[rID][rK]
		if (errAdd != nil || scores[k-1] >= score) && (errNoAdd != nil || scores[k-1] >= noEdgeScore) {
			break
		}
		if errAdd == nil && score > noEdgeScore {
			scores = append(scores, score)
			traces = append(traces, edgeTrace)
		} else if errNoAdd == nil {
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

// Calculates score for given top node v; returns score and best edge.
// k indicates that the edge being added is the k^th edge.
func (dp *DP) scoreV(v *tree.Node, k int, vCycleDP *cycleDPLookups) (bestScore uint, bestCycleTrace *cycleTrace, err error) {
	if k <= 0 {
		panic("scoreV should never be called with zero or negative k value")
	}
	prevK := k - 1
	bestCycleLen := 0
	vCycleDP.solveK(prevK, dp)
	for _, c := range dp.TreeData.Children[v.Id()] {
		if c.Tip() {
			continue
		}
		score, curCycleTrace, err := dp.scoreUDown(v, vCycleDP, prevK)
		if err != nil {
			continue
		}
		cycleLen := dp.cycleLen(curCycleTrace.branch)
		if score > bestScore || bestCycleTrace == nil || (score == bestScore && cycleLen <= bestCycleLen) {
			bestScore = score
			bestCycleTrace = curCycleTrace
			bestCycleLen = cycleLen
		}
	}
	SubtreePostOrder(v, func(u, otherSubtree *tree.Node) {
		score, curCycleTrace, err := dp.scoreUAcross(u, otherSubtree, v, vCycleDP, prevK)
		if err != nil {
			return
		}
		cycleLen := dp.cycleLen(curCycleTrace.branch)
		if score > bestScore || bestCycleTrace == nil || (score == bestScore && cycleLen <= bestCycleLen) {
			bestScore = score
			bestCycleTrace = curCycleTrace
			bestCycleLen = cycleLen
		}
	})
	if bestCycleTrace == nil {
		return maxVal, nil, ErrNoValidSplit
	}
	return bestScore, bestCycleTrace, nil
}

func (dp *DP) cycleLen(br gr.Branch) int {
	if br.Empty() {
		panic("branch was empty!")
	}
	u, w := br.IDs[gr.Ui], br.IDs[gr.Wi]
	v := dp.TreeData.LCA(u, w)
	length := (dp.TreeData.Depths[u] - dp.TreeData.Depths[v]) + (dp.TreeData.Depths[w] - dp.TreeData.Depths[v]) + 1
	if v == u { // we have to account for the edge above v that our new edge is anchored to
		length += 1
	}
	return length
}

func (dp *DP) scoreUDown(v *tree.Node, vCycleScores *cycleDPLookups, prevK int) (bestScore uint, traceback *cycleTrace, err error) {
	SubtreePreOrder(v, func(w *tree.Node) {
		edgeScore := dp.scoreEdge(v, w, v, v)
		wK, wDown, err := BestSplit(vCycleScores.scores[w.Id()], dp.DP[w.Id()], prevK)
		if err != nil { // no valid split, so we don't consider this edge
			return
		}
		wScore, wTrace := vCycleScores.get(w.Id(), wK)
		score := edgeScore + wScore + dp.DP[w.Id()][wDown]
		if score > bestScore || traceback == nil {
			traceback = &cycleTrace{
				pathW:  wTrace,
				wTrace: &dp.Traceback[w.Id()][wDown],
				branch: gr.Branch{IDs: [2]int{v.Id(), w.Id()}},
			}
			bestScore = score
		}
	})
	if traceback == nil {
		return maxVal, nil, ErrNoValidSplit
	}
	return bestScore, traceback, nil
}

// Score branch u -> w (for w in subtree under sub); returns score, best w
func (dp *DP) scoreUAcross(u, sub, v *tree.Node, vCycleScores *cycleDPLookups, prevK int) (bestScore uint, traceback *cycleTrace, err error) {
	if v == u {
		panic("u should not equal v, use scoreUDown instead")
	}
	SubtreePreOrder(sub, func(w *tree.Node) {
		if u == w {
			panic("u should not equal w")
		}
		edgeScore := dp.scoreEdge(u, w, v, sub)
		indices, err := FourWayBestSplit(
			[4][]uint{
				vCycleScores.scores[w.Id()],
				vCycleScores.scores[u.Id()],
				dp.DP[w.Id()],
				dp.DP[u.Id()],
			},
			prevK,
		)
		if err != nil { // no valid split, so we don't consider this edge
			return
		}
		// TODO: fix variable names around here -- what the hell is even this!!
		wK, uK, downW, downU := indices[0], indices[1], indices[2], indices[3]
		wScore, wTrace := vCycleScores.get(w.Id(), wK)
		uScore, uTrace := vCycleScores.get(u.Id(), uK)
		score := edgeScore + wScore + uScore + dp.DP[w.Id()][downW] + dp.DP[u.Id()][downU]
		if score > bestScore || traceback == nil {
			traceback = &cycleTrace{
				pathW:  wTrace,
				pathU:  uTrace,
				wTrace: &dp.Traceback[w.Id()][downW],
				uTrace: &dp.Traceback[u.Id()][downU],
				branch: gr.Branch{IDs: [2]int{u.Id(), w.Id()}},
			}
			bestScore = score
		}
	})
	if traceback == nil {
		return maxVal, nil, ErrNoValidSplit
	}
	return bestScore, traceback, nil
}

func (dp *DP) scoreEdge(u, w, v, wSub *tree.Node) uint {
	score := dp.brScoreCache[u.Id()][w.Id()]
	if score != maxVal {
		return score
	}
	score = uint(0)
	for _, q := range dp.TreeData.Quartets(v.Id()) {
		if QuartetScore(q, u, w, v, wSub, dp.TreeData) == gr.Qeq {
			score += dp.TreeData.NumQuartet(q)
		}
	}
	dp.brScoreCache[u.Id()][w.Id()] = score
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
		switch {
		case !td.InLeafset(v.Id(), t):
			lca = 0
		case td.InLeafset(wSub.Id(), t) || td.InLeafset(u.Id(), bottom):
			lca = td.LCA(w.Id(), tID)
		default:
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

func (dp *DP) traceback(k int) []gr.Branch {
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
