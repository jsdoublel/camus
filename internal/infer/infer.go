package infer

import (
	"fmt"
	"log"

	"github.com/evolbioinfo/gotree/tree"

	gr "github.com/jsdoublel/camus/internal/graphs"
	pr "github.com/jsdoublel/camus/internal/prep"
	sc "github.com/jsdoublel/camus/internal/score"
)

type InferOptions struct {
	NProcs      int                     // number of parallel processes
	QuartetOpts pr.QuartetFilterOptions // quartet filter options
	ScoreMode   sc.ScoreMode            // type of edge score
}

type dpRunner interface {
	RunDP() [][]gr.Branch
}

// Runs Infer algorithm -- returns preprocessed tree data struct, quartet count stats, list of branches.
// Errors returned come from preprocessing (invalid inputs, etc.).
func Infer(tre *tree.Tree, geneTrees []*tree.Tree, opts InferOptions) (*gr.TreeData, [][]gr.Branch, error) {
	log.Println("beginning data preprocessing")
	td, err := pr.Preprocess(tre, geneTrees, opts.NProcs, opts.QuartetOpts)
	if err != nil {
		return nil, nil, fmt.Errorf("preprocess error: %w", err)
	}
	log.Println("calculating edge scores")
	dp, err := makeDP(td, opts.ScoreMode, opts.NProcs)
	if err != nil {
		return nil, nil, err
	}
	log.Println("preprocessing finished, beginning dp algorithm")
	return td, dp.RunDP(), nil
}

func makeDP(td *gr.TreeData, scoreMode sc.ScoreMode, nproces int) (dpRunner, error) {
	n := len(td.Nodes())
	switch scoreMode {
	case sc.MaxScore:
		edgeScores, err := sc.CalculateEdgeScores(sc.MaximizeScorer{}, td, nproces)
		return &DP[uint64]{
			DP:         make([][]uint64, n),
			Traceback:  make([][]trace, n),
			EdgeScores: edgeScores,
			NumNodes:   n,
			Tree:       td,
		}, err
	case sc.NormScore:
		edgeScores, err := sc.CalculateEdgeScores(&sc.NormalizedScorer{}, td, nproces)
		return &DP[float64]{
			DP:         make([][]float64, n),
			Traceback:  make([][]trace, n),
			EdgeScores: edgeScores,
			NumNodes:   n,
			Tree:       td,
		}, err
	case sc.SymScore:
		edgeScores, err := sc.CalculateEdgeScores(&sc.SymDiffScorer{}, td, nproces)
		return &DP[int64]{
			DP:         make([][]int64, n),
			Traceback:  make([][]trace, n),
			EdgeScores: edgeScores,
			NumNodes:   n,
			Tree:       td,
		}, err
	default:
		panic(fmt.Sprintf("invalid score mode (%d)", scoreMode))
	}
}
