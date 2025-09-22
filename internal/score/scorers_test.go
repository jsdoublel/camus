package score

import (
	"errors"
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
	return gr.MakeTreeData(tre, nil)
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
		{name: "positive", alpha: 0.25},
		{name: "large positive", alpha: 1.25},
		{name: "zero", alpha: 0, wantErr: true},
		{name: "negative", alpha: -0.1, wantErr: true},
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
			n := len(td.Nodes())
			if len(scorer.penalties) != n {
				t.Fatalf("penalties length = %d, want %d", len(scorer.penalties), n)
			}
			for i := range scorer.penalties {
				if len(scorer.penalties[i]) != n {
					t.Fatalf("penalties[%d] length = %d, want %d", i, len(scorer.penalties[i]), n)
				}
			}
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
			n := len(td.Nodes())
			if len(scorer.penalties) != n {
				t.Fatalf("penalties length = %d, want %d", len(scorer.penalties), n)
			}
			for i := range scorer.penalties {
				if len(scorer.penalties[i]) != n {
					t.Fatalf("penalties[%d] length = %d, want %d", i, len(scorer.penalties[i]), n)
				}
			}
		})
	}
}
