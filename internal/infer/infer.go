package infer

import (
	"errors"
	"fmt"
	"log"
	"runtime"

	"github.com/evolbioinfo/gotree/tree"

	gr "github.com/jsdoublel/camus/internal/graphs"
	pr "github.com/jsdoublel/camus/internal/prep"
	sc "github.com/jsdoublel/camus/internal/score"
)

var ErrInvalidOption = errors.New("invalid option combination")

type InferOptions struct {
	NProcs      int                     // number of parallel processes
	QuartetOpts pr.QuartetFilterOptions // quartet filter options
	ScoreMode   sc.ScoreMode            // type of edge score
	Alpha       float64
}

// Interface to make DP struct agnostic to generic type when returned
type dpRunner interface {
	RunDP() [][]gr.Branch
}

func MakeInferOptions(nprocs int, quartOpts pr.QuartetFilterOptions, scoreMode sc.ScoreMode, alpha float64) (*InferOptions, error) {
	if quartOpts.QuartetFilterOff() && alpha != 0 {
		return nil, fmt.Errorf("%w: cannot combine non-zero alpha with -q != 0", ErrInvalidOption)
	}
	return &InferOptions{
		NProcs:      setNProcs(nprocs),
		QuartetOpts: quartOpts,
		ScoreMode:   scoreMode,
		Alpha:       alpha,
	}, nil
}

func setNProcs(nprocs int) int {
	maxProcs := runtime.GOMAXPROCS(0)
	switch {
	case nprocs > maxProcs:
		log.Printf("%d is greater than available processes (%d); limit set to %d\n", nprocs, maxProcs, maxProcs)
		return maxProcs
	case nprocs <= 0:
		log.Printf("number of processes not set; defaulting to %d processes\n", maxProcs)
		return maxProcs
	default:
		return nprocs
	}
}

// Runs Infer algorithm -- returns preprocessed tree data struct, quartet count stats, list of branches.
// Errors returned come from preprocessing (invalid inputs, etc.).
func Infer(tre *tree.Tree, geneTrees []*tree.Tree, opts InferOptions) (*gr.TreeData, [][]gr.Branch, error) {
	log.Println("beginning data preprocessing")
	td, err := pr.Preprocess(tre, geneTrees, opts.NProcs, opts.QuartetOpts)
	if err != nil {
		return nil, nil, fmt.Errorf("preprocess error: %w", err)
	}
	var dp dpRunner
	switch opts.ScoreMode {
	case sc.MaxScore:
		dp, err = newDP(sc.MaximizeScorer{}, td, opts.NProcs)
	case sc.NormScore:
		dp, err = newDP(&sc.NormalizedScorer{NGTree: len(geneTrees)}, td, opts.NProcs)
	case sc.SymScore:
		dp, err = newDP(&sc.SymDiffScorer{Alpha: opts.Alpha}, td, opts.NProcs)
	default:
		panic(fmt.Sprintf("invalid score mode (%d)", opts.ScoreMode))
	}
	if err != nil {
		return nil, nil, err
	}
	log.Println("preprocessing finished, beginning dp algorithm")
	return td, dp.RunDP(), nil
}

// Creates DP struct with appropriate score type
func newDP[S sc.Score](scorer sc.Scorer[S], td *gr.TreeData, nproces int) (*DP[S], error) {
	log.Println("calculating edge scores")
	n := len(td.Nodes())
	if edgeScores, err := sc.CalculateEdgeScores(scorer, td, nproces); err != nil {
		return nil, err
	} else {
		return &DP[S]{
			DP:         make([][]S, n),
			Traceback:  make([][]trace, n),
			EdgeScores: edgeScores,
			NumNodes:   n,
			Tree:       td,
		}, nil
	}
}
