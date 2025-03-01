package prep

import (
	"errors"
	"fmt"
	"log"

	"github.com/evolbioinfo/gotree/tree"
)

/*
	preprocess input data

args;

	*tree.Tree   tre         (input constraint tree)
	[]*tree.Tree geneTrees   (input gene trees)

returns:

	[]*Quartet (list of quartets)
	[][]uint   (LCA matrix)
	[][]uint   (leaf sets)
*/
func Preprocess(tre *tree.Tree, geneTrees []*tree.Tree) (*TreeData, error) {
	tre.UpdateTipIndex()
	if !IsBinary(tre.Root()) {
		return nil, errors.New("Constraint tree is not binary")
	}
	quartets, err := processQuartets(geneTrees, tre)
	if err != nil {
		return nil, err
	}
	treeData := PreprocessTreeData(tre, quartets)
	return treeData, nil
}

func processQuartets(geneTrees []*tree.Tree, tre *tree.Tree) ([]*Quartet, error) {
	treeQuartets, err := QuartetsFromTree(tre.Clone(), tre)
	if err != nil {
		panic(err)
	}
	qSet := make(map[Quartet]bool)
	countTotal := len(geneTrees)
	countNew := 0
	countRetained := 0
	for _, gt := range geneTrees {
		gt.UpdateTipIndex()
		newQuartets, err := QuartetsFromTree(gt, tre)
		if err != nil {
			return nil, err
		}
		for q, b := range newQuartets {
			if b && !treeQuartets[q] {
				countNew++
				qSet[q] = true
			}
		}
	}
	quartets := make([]*Quartet, 0, len(qSet)) // convert to slice for better performance
	for q := range qSet {
		countRetained++
		quartets = append(quartets, &q)
	}
	// log.Printf("%d quartets provided, %d were not in constraint tree, and from those there were %d unique topologie(s)\n", countTotal, countNew, countRetained)
	log.Printf("%d quartets provided, %d were not in constraint tree\n", countTotal, countNew)
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
