package score

import (
	"errors"
	"fmt"
	"strings"
	"testing"

	"github.com/evolbioinfo/gotree/io/newick"

	gr "github.com/jsdoublel/camus/internal/graphs"
)

func makeTreeData(t *testing.T, nwk string) *gr.TreeData {
	t.Helper()
	tre, err := newick.NewParser(strings.NewReader(nwk)).Parse()
	if err != nil {
		t.Fatalf("invalid newick in test: %v", err)
	}
	if err := tre.UpdateTipIndex(); err != nil {
		t.Fatalf("failed to update tip index: %v", err)
	}
	tips := tre.AllTipNames()
	qCounts := make(map[gr.Quartet]uint32)
	if len(tips) >= 4 {
		patterns := []string{
			"((%s,%s),(%s,%s));",
			"((%s,%s),(%s,%s));",
		}
		topo := [][4]int{{0, 1, 2, 3}, {0, 2, 1, 3}}
		for i, pattern := range patterns {
			nw := fmt.Sprintf(pattern, tips[topo[i][0]], tips[topo[i][1]], tips[topo[i][2]], tips[topo[i][3]])
			qTree, err := newick.NewParser(strings.NewReader(nw)).Parse()
			if err != nil {
				t.Fatalf("invalid quartet newick in helper: %v", err)
			}
			quartet, err := gr.NewQuartet(qTree, tre)
			if err != nil {
				t.Fatalf("failed to map quartet: %v", err)
			}
			qCounts[quartet] = 1
		}
	}
	return gr.MakeTreeData(tre, qCounts)
}

func TestParseScorerMap(t *testing.T) {
	testCases := []struct {
		name      string
		key       string
		typeCheck func(InitableScorer) bool
	}{
		{
			name: "max",
			key:  "max",
			typeCheck: func(i InitableScorer) bool {
				_, ok := i.(*MaximizeScorer)
				return ok
			},
		},
		{
			name: "norm",
			key:  "norm",
			typeCheck: func(i InitableScorer) bool {
				_, ok := i.(*NormalizedScorer)
				return ok
			},
		},
		{
			name: "sym",
			key:  "sym",
			typeCheck: func(i InitableScorer) bool {
				_, ok := i.(*SymDiffScorer)
				return ok
			},
		},
	}
	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			scorer, ok := ParseScorer[tc.key]
			if !ok {
				t.Fatalf("missing scorer for key %s", tc.key)
			}
			if !tc.typeCheck(scorer) {
				t.Fatalf("scorer for key %s has unexpected type %T", tc.key, scorer)
			}
		})
	}
	if _, ok := ParseScorer["bad"]; ok {
		t.Fatalf("unexpected scorer registered for invalid key")
	}
}

func TestWithNGtrees(t *testing.T) {
	testCases := []struct {
		name    string
		count   int
		wantErr bool
	}{
		{name: "positive", count: 3},
		{name: "zero", count: 0, wantErr: true},
		{name: "negative", count: -2, wantErr: true},
	}
	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			var opts scorerOpts
			err := WithNGtrees(tc.count)(&opts)
			switch {
			case tc.wantErr && err == nil:
				t.Fatalf("expected error for count %d", tc.count)
			case !tc.wantErr && err != nil:
				t.Fatalf("unexpected error: %v", err)
			case !tc.wantErr && opts.nGTrees != tc.count:
				t.Fatalf("nGTrees = %d, want %d", opts.nGTrees, tc.count)
			}
		})
	}
}

func TestWithAlpha(t *testing.T) {
	testCases := []struct {
		name    string
		alpha   float64
		wantErr bool
	}{
		{name: "positive", alpha: 0.1},
		{name: "one", alpha: 1},
		{name: "large positive", alpha: 1.1, wantErr: true},
		{name: "zero", alpha: 0, wantErr: true},
		{name: "negative", alpha: -1, wantErr: true},
	}
	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			var opts scorerOpts
			err := WithAlpha(tc.alpha)(&opts)
			switch {
			case tc.wantErr && err == nil:
				t.Fatalf("expected error for alpha %f", tc.alpha)
			case !tc.wantErr && err != nil:
				t.Fatalf("unexpected error: %v", err)
			case !tc.wantErr && opts.alpha != tc.alpha:
				t.Fatalf("alpha = %f, want %f", opts.alpha, tc.alpha)
			}
		})
	}
}

func TestMaximizeScorerInit(t *testing.T) {
	testCases := []struct {
		name   string
		tree   string
		nprocs int
	}{
		{name: "basic", tree: "((A,B)a,(C,D)b)r;", nprocs: 2},
	}
	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			td := makeTreeData(t, tc.tree)
			scorer := &MaximizeScorer{}
			if err := scorer.Init(td, tc.nprocs); err != nil {
				t.Fatalf("unexpected error: %v", err)
			}
			if !verifyQuartetTotals(t, td, scorer.quartetTotals) {
				t.Fatalf("expected non-zero quartet totals for %s", tc.name)
			}
		})
	}
}

func TestNormalizedScorerInit(t *testing.T) {
	td := makeTreeData(t, "((A,B),(C,D));")
	testCases := []struct {
		name    string
		options []ScoreOptions
		wantErr bool
		wantNG  int
	}{
		{name: "valid", options: []ScoreOptions{WithNGtrees(4)}, wantNG: 4},
		{name: "invalid option", options: []ScoreOptions{WithNGtrees(0)}, wantErr: true},
	}
	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			var scorer NormalizedScorer
			err := scorer.Init(td, 2, tc.options...)
			switch {
			case tc.wantErr && err == nil:
				t.Fatalf("expected error for test %s", tc.name)
			case tc.wantErr && !errors.Is(err, ErrInvalidScorerOption):
				t.Fatalf("unexpected error type: %v", err)
			case tc.wantErr:
				return
			case err != nil:
				t.Fatalf("unexpected error: %v", err)
			}
			if scorer.NGTree != tc.wantNG {
				t.Fatalf("NGTree = %d, want %d", scorer.NGTree, tc.wantNG)
			}
			if !verifyQuartetTotals(t, td, scorer.quartetTotals) {
				t.Fatalf("expected non-zero quartet totals for %s", tc.name)
			}
			verifyPenalties(t, td, scorer.penalties)
		})
	}
}

func TestSymDiffScorerInit(t *testing.T) {
	td := makeTreeData(t, "((A,B),(C,D));")
	testCases := []struct {
		name    string
		options []ScoreOptions
		alpha   float64
		wantErr bool
	}{
		{name: "valid", options: []ScoreOptions{WithAlpha(0.2)}, alpha: 0.2},
		{name: "invalid option", options: []ScoreOptions{WithAlpha(0)}, wantErr: true},
	}
	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			scorer := &SymDiffScorer{}
			err := scorer.Init(td, 2, tc.options...)
			switch {
			case tc.wantErr && err == nil:
				t.Fatalf("expected error for test %s", tc.name)
			case tc.wantErr && !errors.Is(err, ErrInvalidScorerOption):
				t.Fatalf("unexpected error type: %v", err)
			case tc.wantErr:
				return
			case err != nil:
				t.Fatalf("unexpected error: %v", err)
			}
			if scorer.Alpha != tc.alpha {
				t.Fatalf("Alpha = %f, want %f", scorer.Alpha, tc.alpha)
			}
			if !verifyQuartetTotals(t, td, scorer.quartetTotals) {
				t.Fatalf("expected non-zero quartet totals for %s", tc.name)
			}
			verifyPenalties(t, td, scorer.penalties)
		})
	}
}

func verifyQuartetTotals(t *testing.T, td *gr.TreeData, totals [][]uint64) bool {
	t.Helper()
	n := len(td.Nodes())
	if len(totals) != n {
		t.Fatalf("quartet totals length = %d, want %d", len(totals), n)
	}
	nodes := td.Nodes()
	positive := false
	for u := range totals {
		if len(totals[u]) != n {
			t.Fatalf("quartet totals row %d length = %d, want %d", u, len(totals[u]), n)
		}
		for w := range totals[u] {
			val := totals[u][w]
			if ShouldCalcEdge(u, w, td) {
				if val > 0 {
					positive = true
				}
				continue
			}
			if val != 0 {
				t.Fatalf("unexpected quartet total for edge %s->%s", nodes[u].Name(), nodes[w].Name())
			}
		}
	}
	return positive
}

func verifyPenalties(t *testing.T, td *gr.TreeData, penalties [][]uint64) {
	t.Helper()
	n := len(td.Nodes())
	if len(penalties) != n {
		t.Fatalf("penalties length = %d, want %d", len(penalties), n)
	}
	for u := range penalties {
		if len(penalties[u]) != n {
			t.Fatalf("penalties row %d length = %d, want %d", u, len(penalties[u]), n)
		}
		for w := range penalties[u] {
			if ShouldCalcEdge(u, w, td) && penalties[u][w] == 0 {
				t.Fatalf("expected penalty for edge %d->%d", u, w)
			}
		}
	}
}
