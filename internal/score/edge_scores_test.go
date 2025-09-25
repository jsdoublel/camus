package score

import (
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
		nprocs   int
	}{
		{
			name:     "basic",
			tree:     "((A,B)a,(C,D)b)r;",
			quartets: []quartetCount{{nwk: "((A,C),(B,D));", count: 5}},
			nprocs:   3,
		},
		{
			name: "long",
			tree: "(((A,B)a,(C,D)b)e,(E,(F,G)f)c)r;",
			quartets: []quartetCount{
				{nwk: "((A,E),(B,F));", count: 7},
				{nwk: "((A,F),(B,E));", count: 4},
			},
			nprocs: 2,
		},
	}
	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			td := makeTreeDataWithQuartets(t, tc.tree, tc.quartets)
			qt := &QuartetTotals{}
			if err := qt.CalculateQuartetTotals(td, tc.nprocs); err != nil {
				t.Fatalf("unexpected error: %v", err)
			}
			positive := 0
			for u := range qt.quartetTotals {
				for w := range qt.quartetTotals[u] {
					got := qt.quartetTotals[u][w]
					if ShouldCalcEdge(u, w, td) {
						want := quartetsTotal(u, w, td)
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
		uLabel string
		wLabel string
		want   uint64
	}{
		{name: "match", td: td, uLabel: "A", wLabel: "C", want: 5},
		{name: "reverse", td: td, uLabel: "C", wLabel: "A", want: 5},
		{name: "mismatch", td: td, uLabel: "A", wLabel: "B", want: 0},
		{name: "long cycle match", td: tdLong, uLabel: "A", wLabel: "E", want: 7},
		{name: "long cycle reverse", td: tdLong, uLabel: "E", wLabel: "A", want: 7},
		{name: "long cycle alt", td: tdLong, uLabel: "A", wLabel: "F", want: 4},
		{name: "long cycle mismatch", td: tdLong, uLabel: "B", wLabel: "C", want: 0},
	}
	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			uID := nodeIDByLabel(t, tc.td, tc.uLabel)
			wID := nodeIDByLabel(t, tc.td, tc.wLabel)
			got := quartetsTotal(uID, wID, tc.td)
			if got != tc.want {
				t.Fatalf("quartetsTotal(%s,%s) = %d, want %d", tc.uLabel, tc.wLabel, got, tc.want)
			}
		})
	}
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
