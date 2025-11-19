// Package used for preprocessing necessary data for the CAMUS algorithm
package prep

import (
	"context"
	"errors"
	"fmt"
	"log"
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

// Preprocess necessary data. Returns an error if the constraint tree is not valid
// (e.g., not rooted/binary) or if the gene trees are not valid (bad leaf labels).
func Preprocess(tre *tree.Tree, geneTrees []*tree.Tree, nprocs int, opts QuartetFilterOptions, minSupp float64) (*gr.TreeData, error) {
	if err := tre.UpdateTipIndex(); err != nil {
		return nil, fmt.Errorf("constraint tree %w", ErrMulTree)
	}
	if !tre.Rooted() {
		return nil, fmt.Errorf("constraint tree is %w", ErrUnrooted)
	}
	if !TreeIsBinary(tre) {
		return nil, fmt.Errorf("constraint tree is %w", ErrNonBinary)
	}
	if percent := percentNoSupport(geneTrees); percent != 0 && minSupp != 0 {
		log.Printf("WARNING: %.2f%% of gene tree edges do not have support values", percent)
	}
	log.Printf("reading quartets from gene trees")
	qCounts, err := processQuartets(geneTrees, tre, minSupp, nprocs)
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

type quartetShard struct {
	mu     sync.Mutex
	counts map[gr.Quartet]uint32
}

// Returns map containing counts of quartets in input trees (after filtering out
// quartets from constraint tree).
func processQuartets(geneTrees []*tree.Tree, tre *tree.Tree, minSupp float64, nprocs int) (map[gr.Quartet]uint32, error) {
	var missingOnce sync.Once
	const shardBits = 6
	shardCount := 1 << shardBits
	shards := make([]quartetShard, shardCount)
	for i := range shards {
		shards[i].counts = make(map[gr.Quartet]uint32)
	}
	mask := uint64(shardCount - 1)
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
			if b, err := missmatchTaxaSets(gt, tre); err != nil {
				return err
			} else if b {
				missingOnce.Do(func() {
					log.Println("WARNING: missing taxa detected in one or more gene trees;",
						"this may cause issues with some scoring metrics")
				})
			}
			if minSupp != 0 {
				gt.CollapseLowSupport(minSupp, true)
			}
			newQuartets, err := gr.QuartetsFromTree(gt, tre)
			if err != nil {
				return err
			}
			for q, c := range newQuartets {
				shard := &shards[uint64(q)&mask]
				shard.mu.Lock()
				shard.counts[q] += c
				shard.mu.Unlock()
			}
			return nil
		})
	}
	if err := g.Wait(); err != nil {
		return nil, err
	}
	qCounts := make(map[gr.Quartet]uint32)
	for i := range shards {
		for q, c := range shards[i].counts {
			qCounts[q] += c
		}
	}
	return qCounts, nil
}

func missmatchTaxaSets(tre1, tre2 *tree.Tree) (bool, error) {
	n1, err := tre1.NbTips()
	if err != nil {
		return false, err
	}
	n2, err := tre2.NbTips()
	if err != nil {
		return false, err
	}
	return n1 != n2, nil
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

// returns percent of edges without support
func percentNoSupport(trees []*tree.Tree) float64 {
	var total, count int
	for _, t := range trees {
		for _, e := range t.Edges() {
			if e.Right().Tip() {
				continue
			}
			if e.Support() != tree.NIL_SUPPORT {
				count++
			}
			total++
		}
	}
	return float64(total-count) / float64(total) * 100
}
