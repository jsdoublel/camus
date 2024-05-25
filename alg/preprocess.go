package alg

import (
	// "fmt"

	"github.com/evolbioinfo/gotree/tree"
)

/* returns validated quartets and LCA matrix */
func Preprocess(tre *tree.Tree, rawQuartets []*tree.Tree) ([]*tree.Quartet, [][]uint, error) {
	tre.UpdateTipIndex()
	// quartets, err := processQuartets(rawQuartets, tre)
	return nil, nil, nil
}

func processQuartets(rawQuartets []*tree.Tree, tre *tree.Tree) ([]*tree.Quartet, error) {
	quartets := []*tree.Quartet{}
	// for i, rq := range rawQuartets {
	// 	q, err := validateQuartet(rq, tre)
	// 	// check duplicate
	//
	// }
	return quartets, nil
}

/*
checks that newick tree is a valid input quartet and returns *tree.Quartet

- Errors:
  - Does not have four leaves
  - Isn't formatted correctly (does not contain bipartition)
  - Contains taxa not in our tree

- Invalid (i.e., returns nil):
  - Contained in tree

- Panic!!
  - Reads more/less than one quartet despite checks
*/
func validateQuartet(rawQuartet, tre *tree.Tree) (*Quartet, error) {
	// check four leaves
	// TODO: validate all tips are in tre

	// convert to type Quartet
	// qTaxa := make([]uint, 4)
	// count := 0                  // used validate that we only read one quartet from each tree
	// rawQuartet.UpdateTipIndex()
	// rawQuartet.Quartets(true, func(q *tree.Quartet) {
	// 	count++
	// 	qTaxa[0] = q.taxa[0]; qTaxa[1] = q.taxa[1]; qTaxa[2] = q.taxa[2]; qTaxa[3] = q.taxa[3]
	// })
	// if count != 1 {
	// 	panic(fmt.Sprintf("Quartet not read properly (read %d)", count))
	// }
	// map tip indices (qTaxa) to our constraint tree indices
	quartet, err := NewQuartet(rawQuartet, tre)
	if err != nil {
		return nil, err
	}

	// check that quartet is not in tree, but on same taxa set
	return quartet, nil
}

// func convertQuartet(rawQuartet, tre *tree.Tree) (*tree.Quartet, error) {
// 	rawQuartet.UnRoot()
// 	root := rawQuartet.Root() // gets root, unlike UnRoot which modifies the tree
// 	if root.Nneigh() != 3 {
// 		return nil, errors.New("Bad quartet, probably does not contain bipartition")
// 	}
// 	var quartet *Quartet
// 	// the traversal just reads the quartet into stuct with the correct topology
// 	rawQuartet.PostOrder(func(cur, prev *tree.Node, e *tree.Edge) (keep bool) {
// 		if b, err := cur.Parent(); (b == root && cur.Tip()) {
// 			if quartet.taxa[0] == 0 {
// 				ti, err := tre.TipIndex(cur.Name())
// 				qConvertPanic(err)
// 				quartet.taxa[0] = uint(ti)
// 			} else {
// 				ti, err := tre.TipIndex(cur.Name())
// 				qConvertPanic(err)
// 				quartet.taxa[1] = uint(ti)
// 			}
// 		} else if cur.Tip() {
// 			if quartet.taxa[2] == 0 {
// 				ti, err := tre.TipIndex(cur.Name())
// 				qConvertPanic(err)
// 				quartet.taxa[2] = uint(ti)
// 			} else {
// 				ti, err := tre.TipIndex(cur.Name())
// 				qConvertPanic(err)
// 				quartet.taxa[3] = uint(ti)
// 			}
// 		} else if err.Error() == "The node has more than one parent" { // we ignore the error produced when cur = root
// 			panic(fmt.Errorf("convertQuartet: %w", err))
// 		}
// 		return true
// 	})
// 	// assert that quartet is set up correctly
// 	zeroCount := 0 // check that values are set (i.e., at most only one zero value)
// 	return quartet, nil
// }
