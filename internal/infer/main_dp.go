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

	gr "github.com/jsdoublel/camus/internal/graphs"
	sc "github.com/jsdoublel/camus/internal/score"
)

var ErrNoValidSplit = errors.New("no valid split")

// Stores main dp algorithm data
type DP[S sc.Score] struct {
	DP        [][]S        // score for each dp subproblem (DP[v][k])
	Traceback [][]trace    // traceback for each dp subproblem (Traceback[v][k])
	Tree      *gr.TreeData // preprocessed data for our constraint tree
	NumNodes  int          // number of nodes
	Scorer    sc.Scorer[S] // scorer
}

// Stores DP info for lookups corresponding to a given vertex v
type cycleDP[S sc.Score] struct {
	v          *tree.Node
	scores     [][]S               // score for each path (scores[w][k]); unique struct exists for each v
	traceNodes [][]*cycleTraceNode // backtrace for each path (traceNodes[w][k])
}

// ----- Internal Cycle DP Code

// Updates the cycle lookup DP struct for values of k up to prevK
func (cdp *cycleDP[S]) update(prevK int, dp *DP[S]) {
	SubtreePreOrder(cdp.v, func(cur *tree.Node) {
		if prevK == 0 {
			cdp.scores[cur.Id()] = make([]S, 0)
			cdp.traceNodes[cur.Id()] = make([]*cycleTraceNode, 0)
		}
		cdp.grow(cur.Id())
		if len(cdp.scores[cur.Id()])-1 != prevK {
			panic(fmt.Sprintf("wrong size cycle dp tables: len %d, k %d", len(cdp.scores), prevK))
		}
		if cur == cdp.v { // don't want to look at parent of root/v
			return
		}
		p, err := cur.Parent()
		if cdp.v != dp.Tree.Root() && err != nil {
			panic(err)
		} else if p == cdp.v { // if parent is v, then sibling node of cur is also in the cycle
			return
		}
		sibId := dp.Tree.Sibling(cur).Id()
		pScores, pTraces := cdp.scores[p.Id()], cdp.traceNodes[p.Id()]
		pK, sibK, err := BestSplit(pScores, dp.DP[sibId], prevK)
		if err != nil {
			return
		}
		cdp.set(
			cur.Id(),
			prevK,
			pScores[pK]+dp.DP[sibId][sibK],
			cycleTraceNode{p: pTraces[pK], sib: &dp.Traceback[sibId][sibK]},
		)
	})
}

func (cdp *cycleDP[S]) grow(i int) {
	cdp.scores[i] = append(cdp.scores[i], 0)
	cdp.traceNodes[i] = append(cdp.traceNodes[i], nil)
}

func (cdp *cycleDP[S]) set(i, k int, score S, traceNode cycleTraceNode) {
	cdp.scores[i][k] = score
	cdp.traceNodes[i][k] = &traceNode
}

func (cdp *cycleDP[S]) get(i, k int) (S, *cycleTraceNode) {
	return cdp.scores[i][k], cdp.traceNodes[i][k]
}

// ----- Main DP Code

func (dp *DP[S]) RunDP() *DPResults {
	dp.Tree.PostOrder(func(v, prev *tree.Node, e *tree.Edge) (keep bool) {
		if !v.Tip() {
			scores, edgeTrace := dp.solve(v)
			dp.DP[v.Id()] = scores
			dp.Traceback[v.Id()] = edgeTrace
		} else {
			dp.DP[v.Id()] = make([]S, 1)
			dp.Traceback[v.Id()] = make([]trace, 1, dp.NumNodes)
			dp.Traceback[v.Id()][0] = &noCycleTrace{}
		}
		return true
	})
	return dp.collateResults()
}

func (dp *DP[S]) collateResults() *DPResults {
	numOptimal := len(dp.DP[dp.Tree.Root().Id()]) - 1
	log.Printf("%d edges identified\n", numOptimal)
	log.Println("beginning traceback")
	branches := make([][]gr.Branch, numOptimal)
	qStat := make([]float64, 0, numOptimal)
	for k := range numOptimal + 1 {
		if k != 0 {
			finalScore := dp.DP[dp.Tree.Root().Id()][k]
			log.Printf("dp scored %v at root with %d edges\n", finalScore, k)
			branches[k-1] = dp.traceback(k)
			if percent, err := dp.Scorer.PercentQuartetSat(branches[k-1], dp.Tree); err == nil {
				log.Printf("%f percent of quartets satisfied", percent)
				qStat = append(qStat, percent)
			} else {
				log.Printf("error calculating percent quartets satisfied %s, this is a bug! please report!", err.Error())
				qStat = append(qStat, -1)
			}
		}
	}
	log.Println("done.")
	return &DPResults{Tree: dp.Tree, Branches: branches, QSatScore: qStat}
}

// Solve DP problem for vertex v for all k until it stops improving
func (dp *DP[S]) solve(v *tree.Node) ([]S, []trace) {
	lID, rID := dp.Tree.Children[v.Id()][0].Id(), dp.Tree.Children[v.Id()][1].Id()
	scores := make([]S, 1, dp.NumNodes) // choice of capacity is a bit arbitrary
	traces := make([]trace, 1, dp.NumNodes)
	scores[0] = dp.DP[lID][0] + dp.DP[rID][0]
	traces[0] = &noCycleTrace{[2]*trace{&dp.Traceback[lID][0], &dp.Traceback[rID][0]}}
	vCycleDP := cycleDP[S]{
		v:          v,
		scores:     make([][]S, dp.NumNodes),
		traceNodes: make([][]*cycleTraceNode, dp.NumNodes),
	}
	for k := 1; ; k++ {
		var score S
		var backtrace trace
		if noEdgeScore, noEdgeTrace, err := dp.scoreNoAddEdgeK(lID, rID, k); err == nil {
			score, backtrace = noEdgeScore, noEdgeTrace
		}
		if edgeScore, edgeTrace, err := dp.scoreAddEdgeK(v, k, &vCycleDP); err == nil && edgeScore > score {
			score, backtrace = edgeScore, edgeTrace
		}
		if backtrace == nil || scores[k-1] >= score {
			break
		}
		scores = append(scores, score)
		traces = append(traces, backtrace)
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

// Calculate score for vertex v assuming we do not add an edge
func (dp *DP[S]) scoreNoAddEdgeK(lId, rId, k int) (score S, backtrace *noCycleTrace, err error) {
	lK, rK, err := BestSplit(dp.DP[lId], dp.DP[rId], k)
	score = dp.DP[lId][lK] + dp.DP[rId][rK]
	backtrace = &noCycleTrace{prevs: [2]*trace{&dp.Traceback[lId][lK], &dp.Traceback[rId][rK]}}
	return
}

// Calculates score for given top node v assuming an edge is added; returns
// score and best edge. k indicates that the edge being added is the k^th edge.
func (dp *DP[S]) scoreAddEdgeK(v *tree.Node, k int, vCycleDP *cycleDP[S]) (bestScore S, bestCycleTrace *cycleTrace, err error) {
	if k <= 0 {
		panic("should never be called with zero or negative k value")
	}
	prevK := k - 1
	bestCycleLen := 0
	vCycleDP.update(prevK, dp)
	for _, c := range dp.Tree.Children[v.Id()] {
		if c.Tip() {
			continue
		}
		curScore, curCycleTrace, err := dp.scoreEdgesDown(v, vCycleDP, prevK)
		if err != nil {
			continue
		}
		cycleLen := sc.CycleLength(curCycleTrace.branch.IDs[gr.Ui], curCycleTrace.branch.IDs[gr.Wi], dp.Tree)
		if curScore > bestScore || bestCycleTrace == nil || (curScore == bestScore && cycleLen <= bestCycleLen) {
			bestScore = curScore
			bestCycleTrace = curCycleTrace
			bestCycleLen = cycleLen
		}
	}
	SubtreePostOrder(v, func(u, otherSubtree *tree.Node) {
		curScore, curCycleTrace, err := dp.scoreEdgesAcross(u, otherSubtree, v, vCycleDP, prevK)
		if err != nil {
			return
		}
		cycleLen := sc.CycleLength(curCycleTrace.branch.IDs[gr.Ui], curCycleTrace.branch.IDs[gr.Wi], dp.Tree)
		if curScore > bestScore || bestCycleTrace == nil || (curScore == bestScore && cycleLen <= bestCycleLen) {
			bestScore = curScore
			bestCycleTrace = curCycleTrace
			bestCycleLen = cycleLen
		}
	})
	if bestCycleTrace == nil {
		return 0, nil, ErrNoValidSplit
	}
	return bestScore, bestCycleTrace, nil
}

// Scores edges for a branch going from v to all ancestors w
func (dp *DP[S]) scoreEdgesDown(v *tree.Node, vCycleDP *cycleDP[S], prevK int) (bestScore S, traceback *cycleTrace, err error) {
	SubtreePreOrder(v, func(w *tree.Node) {
		if !sc.ShouldCalcEdge(v.Id(), w.Id(), dp.Tree) {
			return
		}
		edgeScore := dp.Scorer.CalcScore(v.Id(), w.Id(), dp.Tree)
		wPathK, wDownK, err := BestSplit(vCycleDP.scores[w.Id()], dp.DP[w.Id()], prevK)
		if err != nil { // no valid split, so we don't consider this edge
			return
		}
		wScore, wPathTrace := vCycleDP.get(w.Id(), wPathK)
		score := edgeScore + wScore + dp.DP[w.Id()][wDownK]
		if score > bestScore || traceback == nil {
			traceback = &cycleTrace{
				pathW:      wPathTrace,
				wDownTrace: &dp.Traceback[w.Id()][wDownK],
				branch:     gr.Branch{IDs: [2]int{v.Id(), w.Id()}},
			}
			bestScore = score
		}
	})
	if traceback == nil {
		return 0, nil, ErrNoValidSplit
	}
	return bestScore, traceback, nil
}

// Score branch u -> w (for all w in subtree under sub)
func (dp *DP[S]) scoreEdgesAcross(u, sub, v *tree.Node, vCycleDP *cycleDP[S], prevK int) (bestScore S, traceback *cycleTrace, err error) {
	if v == u {
		panic("u should not equal v, use scoreUDown instead")
	}
	SubtreePreOrder(sub, func(w *tree.Node) {
		if u == w {
			panic("u should not equal w")
		}
		edgeScore := dp.Scorer.CalcScore(u.Id(), w.Id(), dp.Tree)
		indices, err := FourWayBestSplit(
			[4][]S{
				vCycleDP.scores[w.Id()],
				vCycleDP.scores[u.Id()],
				dp.DP[w.Id()],
				dp.DP[u.Id()],
			},
			prevK,
		)
		if err != nil { // no valid split, so we don't consider this edge
			return
		}
		wPathK, uPathK, wDownK, uDownK := indices[0], indices[1], indices[2], indices[3]
		wScore, wPathTrace := vCycleDP.get(w.Id(), wPathK)
		uScore, uPathTrace := vCycleDP.get(u.Id(), uPathK)
		score := edgeScore + wScore + uScore + dp.DP[w.Id()][wDownK] + dp.DP[u.Id()][uDownK]
		if score > bestScore || traceback == nil {
			traceback = &cycleTrace{
				pathW:      wPathTrace,
				pathU:      uPathTrace,
				wDownTrace: &dp.Traceback[w.Id()][wDownK],
				uDownTrace: &dp.Traceback[u.Id()][uDownK],
				branch:     gr.Branch{IDs: [2]int{u.Id(), w.Id()}},
			}
			bestScore = score
		}
	})
	if traceback == nil {
		return 0, nil, ErrNoValidSplit
	}
	return bestScore, traceback, nil
}

func (dp *DP[S]) traceback(k int) []gr.Branch {
	return dp.Traceback[dp.Tree.Root().Id()][k].traceback()
}
