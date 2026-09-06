package score

import (
	"fmt"
	"math/rand"
	"strings"
	"testing"

	"github.com/evolbioinfo/gotree/io/newick"

	gr "github.com/jsdoublel/camus/internal/graphs"
)

func TestCalculateQuartetTotals(t *testing.T) {
	testCases := []struct {
		name     string
		tree     string
		quartets []quartetCount
		asSet    bool
		nprocs   int
	}{
		{
			name:     "basic",
			tree:     "((A,B)a,(C,D)b)r;",
			quartets: []quartetCount{{nwk: "((A,C),(B,D));", count: 5}},
			asSet:    false,
			nprocs:   3,
		},
		{
			name: "long",
			tree: "(((A,B)a,(C,D)b)e,(E,(F,G)f)c)r;",
			quartets: []quartetCount{
				{nwk: "((A,E),(B,F));", count: 7},
				{nwk: "((A,F),(B,E));", count: 4},
			},
			asSet:  false,
			nprocs: 2,
		},
	}
	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			td := makeTreeDataWithQuartets(t, tc.tree, tc.quartets)
			qt := &QuartetTotals{}
			if err := qt.CalculateQuartetTotals(td, tc.asSet, tc.nprocs); err != nil {
				t.Fatalf("unexpected error: %v", err)
			}
			positive := 0
			for u := range qt.quartetTotals {
				for w := range qt.quartetTotals[u] {
					got := qt.quartetTotals[u][w]
					if ShouldCalcEdge(u, w, td) {
						want := quartetsTotal(u, w, td, tc.asSet)
						if got != want {
							t.Fatalf("quartetTotals[%d][%d] = %d, want %d", u, w, got, want)
						}
						if want > 0 {
							positive++
						}
					} else if got != 0 {
						uName := td.Nodes()[u].Name()
						wName := td.Nodes()[w].Name()
						t.Fatalf("unexpected non-zero total for edge %d->%d (%s->%s): %d", u, w, uName, wName, got)
					}
				}
			}
			if positive == 0 {
				t.Fatalf("expected at least one positive total in test data")
			}
		})
	}
}

func TestShouldCalcEdge(t *testing.T) {
	td := makeTreeData(t, "((A,(B,C)b)a,(D,E)c)r;")
	testCases := []struct {
		name     string
		uLabel   string
		wLabel   string
		expected bool
	}{
		{name: "valid", uLabel: "b", wLabel: "c", expected: true},
		{name: "ancestor", uLabel: "b", wLabel: "a", expected: false},
		{name: "short cycle", uLabel: "a", wLabel: "c", expected: false},
		{name: "root", uLabel: "r", wLabel: "b", expected: false},
	}
	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			uID := nodeIDByLabel(t, td, tc.uLabel)
			wID := nodeIDByLabel(t, td, tc.wLabel)
			got := ShouldCalcEdge(uID, wID, td)
			if got != tc.expected {
				t.Fatalf("shouldCalcEdge(%s,%s) = %t, want %t", tc.uLabel, tc.wLabel, got, tc.expected)
			}
		})
	}
}

func TestCycleLength(t *testing.T) {
	td := makeTreeData(t, "((A,(B,C)b)a,(D,E)c)r;")
	testCases := []struct {
		name   string
		uLabel string
		wLabel string
		length int
	}{
		{name: "b-c", uLabel: "b", wLabel: "c", length: 4},
		{name: "a-c", uLabel: "a", wLabel: "c", length: 3},
		{name: "A-E", uLabel: "A", wLabel: "E", length: 5},
		{name: "b-B", uLabel: "b", wLabel: "B", length: 3},
	}
	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			uID := nodeIDByLabel(t, td, tc.uLabel)
			wID := nodeIDByLabel(t, td, tc.wLabel)
			got := CycleLength(uID, wID, td)
			if got != tc.length {
				t.Fatalf("CycleLength(%s,%s) = %d, want %d", tc.uLabel, tc.wLabel, got, tc.length)
			}
		})
	}
}

func TestQuartetsTotal(t *testing.T) {
	td := makeTreeDataWithQuartets(t, "((A,B)a,(C,D)b)r;", []quartetCount{{nwk: "((A,C),(B,D));", count: 5}})
	tdLong := makeTreeDataWithQuartets(t, "(((A,B)a,(C,D)b)e,(E,(F,G)f)c)r;", []quartetCount{
		{nwk: "((A,E),(B,F));", count: 7},
		{nwk: "((A,F),(B,E));", count: 4},
	})
	testCases := []struct {
		name   string
		td     *gr.TreeData
		asSet  bool
		uLabel string
		wLabel string
		want   uint64
	}{
		{name: "match", td: td, asSet: false, uLabel: "A", wLabel: "C", want: 5},
		{name: "reverse", td: td, asSet: false, uLabel: "C", wLabel: "A", want: 5},
		{name: "mismatch", td: td, asSet: false, uLabel: "A", wLabel: "B", want: 0},
		{name: "long cycle match", td: tdLong, asSet: false, uLabel: "A", wLabel: "E", want: 7},
		{name: "long cycle reverse", td: tdLong, asSet: false, uLabel: "E", wLabel: "A", want: 7},
		{name: "long cycle alt", td: tdLong, asSet: false, uLabel: "A", wLabel: "F", want: 4},
		{name: "long cycle mismatch", td: tdLong, asSet: false, uLabel: "B", wLabel: "C", want: 0},
	}
	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			uID := nodeIDByLabel(t, tc.td, tc.uLabel)
			wID := nodeIDByLabel(t, tc.td, tc.wLabel)
			got := quartetsTotal(uID, wID, tc.td, tc.asSet)
			if got != tc.want {
				t.Fatalf("quartetsTotal(%s,%s) = %d, want %d", tc.uLabel, tc.wLabel, got, tc.want)
			}
		})
	}
}

func makeBalancedNewick(start, end int) string {
	if start == end {
		return fmt.Sprintf("L%d", start)
	}
	mid := start + (end-start)/2
	left := makeBalancedNewick(start, mid)
	right := makeBalancedNewick(mid+1, end)
	return fmt.Sprintf("(%s,%s)", left, right)
}

func BenchmarkQuartetScore(b *testing.B) {
	nwk := makeBalancedNewick(1, 1000) + ";"
	tre, err := newick.NewParser(strings.NewReader(nwk)).Parse()
	if err != nil {
		b.Fatalf("invalid tree: %v", err)
	}
	if err := tre.UpdateTipIndex(); err != nil {
		b.Fatalf("failed to update tip index: %v", err)
	}
	td := gr.MakeTreeData(tre, nil)
	type benchmarkQuery struct {
		q    gr.Quartet
		u    int
		w    int
		v    int
		wSub int
	}
	rng := rand.New(rand.NewSource(42))
	queries := make(map[int][]benchmarkQuery)
	targetCount := 1000
	counts := make(map[int]int)
	nodes := td.Nodes()
	nNodes := len(nodes)
	for counts[gr.Qeq] < targetCount || counts[gr.Qneq] < targetCount || counts[gr.Qdiff] < targetCount {
		u := rng.Intn(nNodes)
		w := rng.Intn(nNodes)
		if !ShouldCalcEdge(u, w, td) {
			continue
		}
		leaves := [4]uint16{
			uint16(rng.Intn(1000)),
			uint16(rng.Intn(1000)),
			uint16(rng.Intn(1000)),
			uint16(rng.Intn(1000)),
		}
		if hasDuplicates(leaves) {
			continue
		}
		var qTaxa [4]int16
		for i := 0; i < 4; i++ {
			qTaxa[i] = int16(leaves[i])
		}
		v := td.LCA(u, w)
		wSub := getWSubtree(u, w, v, td)
		uNode := td.IdToNodes[u]
		wNode := td.IdToNodes[w]
		vNode := td.IdToNodes[v]
		for _, topo := range []uint8{gr.Qtopo1, gr.Qtopo2, gr.Qtopo3} {
			q := makeBenchmarkQuartet(qTaxa, topo)
			res := QuartetScore(q, uNode, wNode, vNode, wSub, td)
			if counts[res] < targetCount {
				queries[res] = append(queries[res], benchmarkQuery{
					q:    q,
					u:    u,
					w:    w,
					v:    v,
					wSub: wSub.Id(),
				})
				counts[res]++
			}
		}
	}
	testCases := []struct {
		name string
		code int
	}{
		{name: "eq", code: gr.Qeq},
		{name: "neq", code: gr.Qneq},
		{name: "diff", code: gr.Qdiff},
	}
	for _, tc := range testCases {
		b.Run(tc.name, func(b *testing.B) {
			qList := queries[tc.code]
			var got int
			b.ReportAllocs()
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				query := qList[i%len(qList)]
				got = QuartetScore(query.q, td.IdToNodes[query.u], td.IdToNodes[query.w], td.IdToNodes[query.v], td.IdToNodes[query.wSub], td)
			}
			b.StopTimer()
			if got != tc.code {
				b.Fatalf("QuartetScore returned unexpected result: %d", got)
			}
		})
	}
}

func hasDuplicates(leaves [4]uint16) bool {
	return leaves[0] == leaves[1] || leaves[0] == leaves[2] || leaves[0] == leaves[3] ||
		leaves[1] == leaves[2] || leaves[1] == leaves[3] ||
		leaves[2] == leaves[3]
}

func makeBenchmarkQuartet(taxa [4]int16, topo uint8) gr.Quartet {
	var q uint64
	for i, t := range taxa {
		q |= uint64(t) << (15 * i)
	}
	q |= uint64(topo) << 60
	return gr.Quartet(q)
}

type quartetCount struct {
	nwk   string
	count uint32
}

func makeTreeDataWithQuartets(t *testing.T, treeNWK string, quartets []quartetCount) *gr.TreeData {
	t.Helper()
	tre, err := newick.NewParser(strings.NewReader(treeNWK)).Parse()
	if err != nil {
		t.Fatalf("invalid tree newick: %v", err)
	}
	if err := tre.UpdateTipIndex(); err != nil {
		t.Fatalf("failed to update tip index: %v", err)
	}
	qCounts := make(map[gr.Quartet]uint32)
	for _, qt := range quartets {
		qTree, err := newick.NewParser(strings.NewReader(qt.nwk)).Parse()
		if err != nil {
			t.Fatalf("invalid quartet newick %s: %v", qt.nwk, err)
		}
		q, err := gr.NewQuartet(qTree, tre)
		if err != nil {
			t.Fatalf("failed to build quartet %s: %v", qt.nwk, err)
		}
		qCounts[q] = qt.count
	}
	return gr.MakeTreeData(tre, qCounts)
}
