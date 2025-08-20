// Package used for preprocessing necessary data for the CAMUS algorithm
package prep

import (
	"cmp"
	"context"
	"errors"
	"fmt"
	"log"
	"slices"
	"strconv"
	"sync"

	"github.com/evolbioinfo/gotree/tree"
	"golang.org/x/sync/errgroup"

	gr "github.com/jsdoublel/camus/internal/graphs"
)

var (
	ErrUnrooted     = errors.New("not rooted")
	ErrNonBinary    = errors.New("not binary")
	ErrMulTree      = errors.New("contains duplicate labels")
	ErrTypeOutRange = errors.New("out of type range")
)

// Options for quartet filter mode
type QuartetFilterOptions struct {
	mode      QMode     // mode (value between 0 and 1)
	threshold Threshold // threshold for filtering [0, 1]
}

func SetQuartetFilterOptions(mode int, threshold float64) (*QuartetFilterOptions, error) {
	var m QMode
	if err := m.Set(mode); err != nil {
		return nil, err
	}
	var t Threshold
	if err := t.Set(threshold); err != nil {
		return nil, err
	}
	return &QuartetFilterOptions{mode: m, threshold: t}, nil
}

func (opts QuartetFilterOptions) QuartetFilterOff() bool {
	return opts.mode == 0
}

type QMode int

const (
	NonRestrictive QMode = iota + 1
	Restrictive
)

func (mode *QMode) Set(n int) error {
	if n < 0 || n > 2 {
		return fmt.Errorf("quartet mode %d is %w", n, ErrTypeOutRange)
	}
	*mode = QMode(n)
	return nil
}

func (mode QMode) String() string {
	return strconv.Itoa(int(mode))
}

type Threshold float64

func (thresh *Threshold) Set(n float64) error {
	if n < 0 || n > 1 {
		return fmt.Errorf("threshold %f is %w", n, ErrTypeOutRange)
	}
	*thresh = Threshold(n)
	return nil
}

func (thresh Threshold) String() string {
	return strconv.FormatFloat(float64(thresh), 'f', -1, 64)
}

func (thresh Threshold) Keep(counts []uint32) bool {
	if len(counts) != 3 {
		panic("there should be three counts, one for each quartet topology")
	}
	slices.Sort(counts)
	sum := counts[0] + counts[1]
	return uint32(float64(thresh)*float64(sum)) < counts[1]-counts[0]
}

func filterQuartets(qCounts map[gr.Quartet]uint32, opts QuartetFilterOptions) {
	for q := range qCounts {
		quartets := q.AllQuartets()
		counts := []uint32{qCounts[quartets[0]], qCounts[quartets[1]], qCounts[quartets[2]]}
		slices.SortFunc(quartets, func(q1, q2 gr.Quartet) int {
			return cmp.Compare(qCounts[q1], qCounts[q2])
		})
		if !opts.threshold.Keep(counts) {
			delete(qCounts, quartets[0])
			delete(qCounts, quartets[1])
			continue
		}
		switch opts.mode {
		case NonRestrictive:
		case Restrictive:
			delete(qCounts, quartets[0])
		default:
			panic("invalid quartet mode case")
		}
	}
}

// Preprocess necessary data. Returns an error if the constraint tree is not valid
// (e.g., not rooted/binary) or if the gene trees are not valid (bad leaf labels).
func Preprocess(tre *tree.Tree, geneTrees []*tree.Tree, nprocs int, opts QuartetFilterOptions) (*gr.TreeData, error) {
	if err := tre.UpdateTipIndex(); err != nil {
		return nil, fmt.Errorf("constraint tree %w", ErrMulTree)
	}
	if !tre.Rooted() {
		return nil, fmt.Errorf("constraint tree is %w", ErrUnrooted)
	}
	if !TreeIsBinary(tre) {
		return nil, fmt.Errorf("constraint tree is %w", ErrNonBinary)
	}
	qCounts, err := processQuartets(geneTrees, tre, nprocs)
	if err != nil {
		return nil, err
	}
	if opts.mode != 0 {
		filterQuartets(qCounts, opts)
	}
	treeQuartets, err := gr.QuartetsFromTree(tre.Clone(), tre)
	if err != nil {
		return nil, err
	}
	for q := range treeQuartets {
		delete(qCounts, q)
	}
	log.Printf("%d gene trees provided, containing %d quartets not in the constraint tree\n", len(geneTrees), len(qCounts))
	log.Printf("analyzing constraint tree")
	treeData := gr.MakeTreeData(tre, qCounts)
	return treeData, nil
}

// Returns map containing counts of quartets in input trees (after filtering out
// quartets from constraint tree).
func processQuartets(geneTrees []*tree.Tree, tre *tree.Tree, nprocs int) (map[gr.Quartet]uint32, error) {
	log.Printf("processing quartets")
	qCounts := make(map[gr.Quartet]uint32)
	var mu sync.Mutex
	g, ctx := errgroup.WithContext(context.Background())
	g.SetLimit(nprocs)
	for i, gt := range geneTrees {
		g.Go(func() error {
			if err := ctx.Err(); err != nil {
				return err
			}
			if err := gt.UpdateTipIndex(); err != nil {
				return fmt.Errorf("gene tree on line %d : %w", i+1, ErrMulTree)
			}
			newQuartets, err := gr.QuartetsFromTree(gt, tre)
			if err != nil {
				return err
			}
			localCounts := make(map[gr.Quartet]uint32)
			for quartet, count := range newQuartets {
				localCounts[quartet] += count
			}
			mu.Lock()
			for q, c := range localCounts {
				qCounts[q] += c
			}
			mu.Unlock()
			return nil
		})
	}
	if err := g.Wait(); err != nil {
		return nil, err
	}
	return qCounts, nil
}

func NetworkIsBinary(ntw *tree.Tree) bool {
	if !ntw.Rooted() {
		return false
	}
	neighbors := ntw.Root().Neigh()
	if len(neighbors) != 2 {
		panic("tree is not rooted (even though it is??)")
	}
	return isBinary(neighbors[0], true) && isBinary(neighbors[1], true)
}

func TreeIsBinary(tre *tree.Tree) bool {
	if !tre.Rooted() {
		return false
	}
	neighbors := tre.Root().Neigh()
	if len(neighbors) != 2 {
		panic("tree is not rooted (even though it is??)")
	}
	return isBinary(neighbors[0], false) && isBinary(neighbors[1], false)
}

func isBinary(node *tree.Node, allowUnifurcations bool) bool {
	if node.Tip() {
		return true
	}
	children := gr.GetChildren(node)
	if node.Nneigh() == 2 && allowUnifurcations {
		return isBinary(children[0], allowUnifurcations)
	}
	if node.Nneigh() != 3 {
		return false
	}
	return isBinary(children[0], allowUnifurcations) && isBinary(children[1], allowUnifurcations)
}
