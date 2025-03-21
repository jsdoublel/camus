package prep

import (
	"errors"
	"fmt"
	"log"

	"github.com/evolbioinfo/gotree/tree"
)

var ErrInvalidTree = errors.New("invalid tree")

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
	if !tre.Rooted() {
		return nil, fmt.Errorf("%w, constraint tree is not rooted", ErrInvalidTree)
	}
	if !TreeIsBinary(tre) {
		return nil, fmt.Errorf("%w, constraint tree is not binary", ErrInvalidTree)
	}
	if !IsSingleCopy(tre) {
		return nil, fmt.Errorf("%w, constraint tree has duplicate labels", ErrInvalidTree)
	}
	qCounts, err := processQuartets(geneTrees, tre)
	if err != nil {
		return nil, err
	}
	treeData := MakeTreeData(tre, qCounts)
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
	for i, gt := range geneTrees {
		if !IsSingleCopy(gt) {
			return nil, fmt.Errorf("%w, gene tree on line %d has duplicate labels", ErrInvalidTree, i)
		}
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
	return len(labels) == len(tre.Tips())
}
