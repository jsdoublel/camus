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
	qCounts, err := processQuartets(geneTrees, tre)
	if err != nil {
		return nil, err
	}
	treeData := PreprocessTreeData(tre, qCounts)
	return treeData, nil
}

func processQuartets(geneTrees []*tree.Tree, tre *tree.Tree) (map[Quartet]uint, error) {
	treeQuartets, err := QuartetsFromTree(tre.Clone(), tre)
	if err != nil {
		panic(err)
	}
	qCounts := make(map[Quartet]uint)
	countTotal := len(geneTrees)
	countNew := 0
	for _, gt := range geneTrees {
		gt.UpdateTipIndex()
		newQuartets, err := QuartetsFromTree(gt, tre)
		if err != nil {
			return nil, err
		}
		for quartet, count := range newQuartets {
			if count < 0 {
				panic(fmt.Sprintf("negative quartet count %d", count))
			}
			if treeQuartets[quartet] == 0 {
				if qCounts[quartet] == 0 {
					countNew++
				}
				qCounts[quartet] += count
			}
		}
	}
	log.Printf("%d gene trees provided, %d new quartet trees were found\n", countTotal, countNew)
	return qCounts, nil
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
