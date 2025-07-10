// Package used for preprocessing necessary data for the CAMUS algorithm
package prep

import (
	"errors"
	"fmt"
	"log"
	"math"
	"slices"
	"strconv"

	"github.com/evolbioinfo/gotree/tree"

	gr "github.com/jsdoublel/camus/internal/graphs"
)

var (
	ErrUnrooted  = errors.New("not rooted")
	ErrNonBinary = errors.New("not binary")
	ErrMulTree   = errors.New("contains duplicate labels")
)

type QMode int

func (mode *QMode) Set(n string) error {
	res, err := strconv.Atoi(n)
	if err != nil {
		return err
	}
	if res < 0 || res > 2 {
		return fmt.Errorf("quartet mode %d out of range", res)
	}
	*mode = QMode(res)
	return nil
}

func (mode *QMode) String() string {
	return strconv.Itoa(int(*mode))
}

// Preprocess necessary data. Returns an error if the constraint tree is not valid
// (e.g., not rooted/binary) or if the gene trees are not valid (bad leaf labels).
func Preprocess(tre *tree.Tree, geneTrees []*tree.Tree, qMode QMode) (*gr.TreeData, error) {
	if err := tre.UpdateTipIndex(); err != nil {
		return nil, fmt.Errorf("constraint tree %w", ErrMulTree)
	}
	if !tre.Rooted() {
		return nil, fmt.Errorf("constraint tree is %w", ErrUnrooted)
	}
	if !TreeIsBinary(tre) {
		return nil, fmt.Errorf("constraint tree is %w", ErrNonBinary)
	}
	qCounts, err := processQuartets(geneTrees, tre, qMode)
	if err != nil {
		return nil, err
	}
	treeData := gr.MakeTreeData(tre, qCounts)
	return treeData, nil
}

// Returns map containing counts of quartets in input trees (after filtering out
// quartets from constraint tree).
func processQuartets(geneTrees []*tree.Tree, tre *tree.Tree, qMode QMode) (map[gr.Quartet]uint, error) {
	treeQuartets, err := gr.QuartetsFromTree(tre.Clone(), tre)
	if err != nil {
		panic(err)
	}
	taxaSets := make(map[[4]int]struct{})
	qCounts := make(map[gr.Quartet]uint)
	countGTree := len(geneTrees)
	countTotal := uint(0)
	countNew := uint(0)
	for i, gt := range geneTrees {
		LogEveryNPercent(i, 10, len(geneTrees), fmt.Sprintf("processed %d out of %d gene trees", i+1, countGTree))
		if err := gt.UpdateTipIndex(); err != nil {
			return nil, fmt.Errorf("gene tree on line %d : %w", i, ErrMulTree)
		}
		newQuartets, err := gr.QuartetsFromTree(gt, tre)
		if err != nil {
			return nil, err
		}
		for quartet, count := range newQuartets {
			if _, exists := taxaSets[quartet.Taxa]; !exists {
				taxaSets[quartet.Taxa] = struct{}{}
			}
			if treeQuartets[quartet] == 0 {
				qCounts[quartet] += count
				countNew += count
			}
			countTotal += count
		}
	}
	log.Printf("%d gene trees provided, containing %d quartets; %d new quartet trees were found\n", countGTree, countTotal, countNew)
	if qMode != 0 {
		filterQuartets(qCounts, taxaSets, qMode)
	}
	return qCounts, nil
}

func filterQuartets(qCounts map[gr.Quartet]uint, taxaSets map[[4]int]struct{}, qMode QMode) {
	for taxaSet := range taxaSets {
		quartets := []gr.Quartet{gr.Quartet{Taxa: taxaSet, Topology: gr.Qtopo1},
			gr.Quartet{Taxa: taxaSet, Topology: gr.Qtopo2},
			gr.Quartet{Taxa: taxaSet, Topology: gr.Qtopo3}}
		counts := []uint{qCounts[quartets[0]], qCounts[quartets[1]], qCounts[quartets[2]]}
		switch qMode {
		case 1:
			minQ := quartets[slices.Index(counts, slices.Min(counts))]
			delete(qCounts, minQ)
		case 2:
			maxQ := quartets[slices.Index(counts, slices.Max(counts))]
			for _, q := range quartets {
				if q != maxQ {
					delete(qCounts, q)
				}
			}
		default:
			panic("invalid quartet mode case")
		}
	}
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

func LogEveryNPercent(i, n, total int, message string) {
	if (i+1)%max(total/int(math.Ceil(float64(100)/float64(n))), 1) == 0 {
		log.Print(message)
	}
}
