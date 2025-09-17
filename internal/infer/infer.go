package infer

import (
	"fmt"
	"log"

	"github.com/evolbioinfo/gotree/tree"

	gr "github.com/jsdoublel/camus/internal/graphs"
	pr "github.com/jsdoublel/camus/internal/prep"
	"github.com/jsdoublel/camus/internal/score"
)

type InferOptions struct {
	NProcs      int                     // number of parallel processes
	QuartetOpts pr.QuartetFilterOptions // quartet filter options
	ScoreType   int                     // type of edge score
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
	// edgeScores, err := score.CalculateEdgeScores(opts.ScoreType, td, opts.NProcs)
	// if err != nil {
	// 	return nil, nil, err
	// }
	dp, err := makeDP(td, opts.ScoreType, opts.NProcs)
	if err != nil {
		return nil, nil, err
	}
	log.Println("preprocessing finished, beginning dp algorithm")
	// dp := &DP{
	// 	DP:         make([][]uint, n),
	// 	Traceback:  make([][]trace, n),
	// 	Tree:       td,
	// 	EdgeScores: edgeScores,
	// 	NumNodes:   n,
	// }
	return td, dp.RunDP(), nil
}

func makeDP(td *gr.TreeData, scoreType, nproces int) (*DP, error) {
	edgeScores, err := score.CalculateEdgeScores(score.MaximizeScorer{}, td, nproces)
	if err != nil {
		return nil, err
	}
	n := len(td.Nodes())
	dp := &DP{
		DP:         make([][]uint64, n),
		Traceback:  make([][]trace, n),
		Tree:       td,
		EdgeScores: edgeScores,
		NumNodes:   n,
	}
	return dp, nil
}
