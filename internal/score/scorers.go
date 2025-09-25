package score

import (
	"errors"
	"fmt"

	gr "github.com/jsdoublel/camus/internal/graphs"
)

var ErrInvalidScorerOption = errors.New("invalid scorer option")

var ParseScorer = map[string]InitableScorer{
	"max":  &MaximizeScorer{},
	"norm": &NormalizedScorer{},
	"sym":  &SymDiffScorer{},
}

// interface to allow scorers to be stored in a map together
type InitableScorer interface {
	Init(td *gr.TreeData, nprocs int, opts ...ScoreOptions) error
}

type ScoreOptions func(opts *scorerOpts) error

type scorerOpts struct {
	nGTrees int
	alpha   int64
}

type Score interface{ int64 | uint64 | float64 }

// scorers implement different scorring metrics
type Scorer[S Score] interface {
	Init(td *gr.TreeData, nprocs int, opts ...ScoreOptions) error
	CalcScore(u, w int, td *gr.TreeData) S
	TotalSatQuartets(branches []gr.Branch) (uint64, error)
}

type MaximizeScorer struct {
	QuartetTotals
}

func (s *MaximizeScorer) Init(td *gr.TreeData, nprocs int, opts ...ScoreOptions) error {
	// No opts are needed for MaximizeScorer
	return s.CalculateQuartetTotals(td, nprocs)
}

func (s MaximizeScorer) CalcScore(u, w int, td *gr.TreeData) uint64 {
	return s.quartetTotals[u][w]
}

type NormalizedScorer struct {
	QuartetTotals
	NGTree    int
	penalties [][]uint64
}

func WithNGtrees(ngtrees int) ScoreOptions {
	return func(options *scorerOpts) error {
		if ngtrees <= 0 {
			return fmt.Errorf("%w, number of gene trees must be positive, but is %d", ErrInvalidScorerOption, ngtrees)
		}
		options.nGTrees = ngtrees
		return nil
	}
}

func (s *NormalizedScorer) Init(td *gr.TreeData, nprocs int, opts ...ScoreOptions) error {
	var options scorerOpts
	for _, opt := range opts {
		if err := opt(&options); err != nil {
			return err
		}
	}
	s.NGTree = options.nGTrees
	if err := s.CalculateQuartetTotals(td, nprocs); err != nil {
		return err
	}
	var err error
	if s.penalties, err = CalculateEdgePenalties(td, nprocs); err != nil {
		return err
	}
	return nil
}

func (s NormalizedScorer) CalcScore(u, w int, td *gr.TreeData) float64 {
	return float64(s.quartetTotals[u][w]) / (float64(s.NGTree) * float64(s.penalties[u][w]))
}

type SymDiffScorer struct {
	QuartetTotals
	Alpha     int64
	penalties [][]uint64
}

func WithAlpha(alpha int64) ScoreOptions {
	return func(options *scorerOpts) error {
		if alpha <= 0 {
			return fmt.Errorf("%w, alpha must be in greater than zero, but is %d", ErrInvalidScorerOption, alpha)
		}
		options.alpha = alpha
		return nil
	}
}

func (s *SymDiffScorer) Init(td *gr.TreeData, nprocs int, opts ...ScoreOptions) error {
	var options scorerOpts
	for _, opt := range opts {
		if err := opt(&options); err != nil {
			return err
		}
	}
	s.Alpha = options.alpha
	if err := s.CalculateQuartetTotals(td, nprocs); err != nil {
		return err
	}
	var err error
	if s.penalties, err = CalculateEdgePenalties(td, nprocs); err != nil {
		return err
	}
	return nil
}

func (s SymDiffScorer) CalcScore(u, w int, td *gr.TreeData) int64 {
	return int64(s.quartetTotals[u][w]) - s.Alpha*int64(s.penalties[u][w])
}
