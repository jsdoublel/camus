package prep

import (
	"errors"
	"fmt"
	"os"

	"github.com/evolbioinfo/gotree/tree"
)

/*
	preprocess input data

args;

	*tree.Tree   tre         (input constraint tree)
	[]*tree.Tree rawQuartets (unprocess quartet trees)

returns:

	[]*Quartet (list of quartets)
	[][]uint   (LCA matrix)
	[][]uint   (leaf sets)
*/
func Preprocess(tre *tree.Tree, rawQuartets []*tree.Tree) ([]*Quartet, *TreeData, error) {
	tre.UpdateTipIndex()
	if !IsBinary(tre) {
		return nil, nil, errors.New("Constraint tree is not binary")
	}
	quartets, err := processQuartets(rawQuartets, tre)
	if err != nil {
		return nil, nil, err
	}
	treeData := PreprocessTreeData(tre)
	return quartets, treeData, nil
}

func processQuartets(rawQuartets []*tree.Tree, tre *tree.Tree) ([]*Quartet, error) {
	treeQuartets := QuartetsFromTree(tre)
	qSet := make(map[Quartet]bool)
	countTotal := len(rawQuartets)
	countNew := 0
	countRetained := 0
	for i, rq := range rawQuartets {
		q, err := validateQuartet(rq, tre)
		if err != nil {
			return nil, fmt.Errorf("Bad quartet line %d: %w", i, err)
		} else if !treeQuartets[*q] {
			countNew++
			qSet[*q] = true
		}
	}
	quartets := make([]*Quartet, 0, len(qSet)) // convert to slice for better performance
	for q := range qSet {
		countRetained++
		quartets = append(quartets, &q)
	}
	fmt.Fprintf(os.Stderr, "%d quartets provided\n%d were not in constraint tree\n%d retained after removing duplicates\n", countTotal, countNew, countRetained)
	return quartets, nil
}

/*
checks that newick tree is a valid input quartet and returns *Quartet

errors:

	Does not have four leaves
	Isn't formatted correctly (does not contain bipartition)
	Contains taxa not in our tree
*/
func validateQuartet(rawQuartet, tre *tree.Tree) (*Quartet, error) {
	for _, leaf := range rawQuartet.AllTipNames() {
		if b, err := tre.ExistsTip(leaf); err != nil {
			panic(fmt.Errorf("validate quartets %w", err))
		} else if !b {
			return nil, errors.New("contains taxa not in constraint tree")
		}
	}
	quartet, err := NewQuartet(rawQuartet, tre)
	if err != nil {
		return nil, err
	}
	return quartet, nil
}
