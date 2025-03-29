// Package used for preprocessing necessary data for the CAMUS algorithm
package prep

import (
	"errors"
	"fmt"
	"log"

	"github.com/evolbioinfo/gotree/tree"

	"github.com/jsdoublel/camus/qrt"
)

var ErrInvalidTree = errors.New("invalid tree")

// Preprocess necessary data. Returns an error if the constraint tree is not valid
// (e.g., not rooted/binary) or if the gene trees are not valid (bad leaf labels).
func Preprocess(tre *tree.Tree, geneTrees []*tree.Tree) (*TreeData, *qrt.QuartetStats, error) {
	tre.UpdateTipIndex()
	if !tre.Rooted() {
		return nil, nil, fmt.Errorf("%w, constraint tree is not rooted", ErrInvalidTree)
	}
	if !TreeIsBinary(tre) {
		return nil, nil, fmt.Errorf("%w, constraint tree is not binary", ErrInvalidTree)
	}
	if !IsSingleCopy(tre) {
		return nil, nil, fmt.Errorf("%w, constraint tree has duplicate labels", ErrInvalidTree)
	}
	qCounts, qStats, err := processQuartets(geneTrees, tre)
	if err != nil {
		return nil, nil, err
	}
	treeData := MakeTreeData(tre, qCounts)
	return treeData, qStats, nil
}

// Returns map containing counts of quartets in input trees (after filtering out
// quartets from constraint tree).
func processQuartets(geneTrees []*tree.Tree, tre *tree.Tree) (*map[qrt.Quartet]uint, *qrt.QuartetStats, error) {
	qStats := qrt.MakeQStats()
	treeQuartets, err := qrt.QuartetsFromTree(tre.Clone(), tre)
	if err != nil {
		panic(err)
	}
	qCounts := make(map[qrt.Quartet]uint)
	countTotal := len(geneTrees)
	countNew := uint(0)
	for i, gt := range geneTrees {
		if !IsSingleCopy(gt) {
			return nil, nil, fmt.Errorf("%w, gene tree on line %d has duplicate labels", ErrInvalidTree, i)
		}
		gt.UpdateTipIndex()
		newQuartets, err := qrt.QuartetsFromTree(gt, tre)
		if err != nil {
			return nil, nil, err
		}
		for quartet, count := range newQuartets {
			if count < 0 {
				panic(fmt.Sprintf("negative quartet count %d", count))
			}
			qStats.Count(&quartet, count)
			if treeQuartets[quartet] == 0 {
				qCounts[quartet] += count
				countNew += count
			}
		}
	}
	log.Printf("%d gene trees provided, %d new quartet trees were found\n", countTotal, countNew)
	return &qCounts, qStats, nil
}

func TreeIsBinary(tre *tree.Tree) bool {
	if !tre.Rooted() {
		return false
	}
	neighbors := tre.Root().Neigh()
	if len(neighbors) != 2 {
		panic("tree is not rooted (even though it is??)")
	}
	return isBinary(neighbors[0]) && isBinary(neighbors[1])
}

func isBinary(node *tree.Node) bool {
	if node.Tip() {
		return true
	}
	if node.Nneigh() != 3 {
		return false
	}
	children := getChildren(node)
	return isBinary(children[0]) && isBinary(children[1])
}

func IsSingleCopy(tre *tree.Tree) bool {
	labels := make(map[string]bool)
	for _, l := range tre.Tips() {
		labels[l.Name()] = true
	}
	// have to use less efficient tre.Tips() because of weird behavior of gotree with multrees
	return len(labels) == len(tre.Tips())
}
