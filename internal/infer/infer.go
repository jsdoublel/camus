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
	ScoreMode   sc.InitableScorer       // type of edge score
	Alpha       int64
}

// Interface to make DP struct agnostic to generic type when returned
type dpRunner interface {
	RunDP() [][]gr.Branch
}

func MakeInferOptions(nprocs int, quartOpts pr.QuartetFilterOptions, scoreMode sc.InitableScorer, alpha int64) (*InferOptions, error) {
	if _, ok := scoreMode.(*sc.SymDiffScorer); !ok && alpha != 0 {
		return nil, fmt.Errorf("%w: cannot combine non-zero alpha with ", ErrInvalidOption)
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
	switch scorer := opts.ScoreMode.(type) {
	case *sc.MaximizeScorer:
		dp, err = newDP(scorer, td, opts.NProcs)
	case *sc.NormalizedScorer:
		dp, err = newDP(scorer, td, opts.NProcs, sc.WithNGtrees(len(geneTrees)))
	case *sc.SymDiffScorer:
		dp, err = newDP(scorer, td, opts.NProcs, sc.WithAlpha(opts.Alpha))
	default:
		panic(fmt.Sprintf("unsupported scorer type %T", scorer))
	}
	if err != nil {
		return nil, nil, err
	}
	log.Println("preprocessing finished, beginning dp algorithm")
	return td, dp.RunDP(), nil
}

// Creates DP struct with appropriate score type
func newDP[S sc.Score](scorer sc.Scorer[S], td *gr.TreeData, nprocs int, opts ...sc.ScoreOptions) (*DP[S], error) {
	if err := scorer.Init(td, nprocs, opts...); err != nil {
		return nil, err
	}
	log.Println("calculating edge scores")
	n := len(td.Nodes())
	return &DP[S]{
		DP:        make([][]S, n),
		Traceback: make([][]trace, n),
		// EdgeScores: edgeScores,
		Scorer:   scorer,
		NumNodes: n,
		Tree:     td,
	}, nil
}
