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
	n := len(td.Nodes())
	edgePenalties, err := score.CalcuateEdgePenalties(td, opts.NProcs)
	if err != nil {
		return nil, nil, err
	}
	edgeScores, err := score.CalculateEdgeScores(td, opts.NProcs)
	if err != nil {
		return nil, nil, err
	}
	log.Println("preprocessing finished, beginning dp algorithm")
	dp := &DP{
		DP:            make([][]uint, n),
		Traceback:     make([][]trace, n),
		Tree:          td,
		EdgeScores:    edgeScores,
		EdgePenalties: edgePenalties,
		NumNodes:      n,
	}
	return td, dp.RunDP(), nil
}
