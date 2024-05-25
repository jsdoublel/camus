package alg

import (
	"errors"
	"fmt"
	"slices"

	"github.com/evolbioinfo/gotree/tree"
)

type Quartet struct {
	taxa     []int // should be in sorted order
	topology uint8
}

/* enumerate possible topologies for a quartet */
const (
	TOPONIL = 0b0000
	TOPO1   = 0b0011
	TOPO2   = 0b0101
	TOPO3   = 0b0110
)

/* possible results for quartet comparison */
const (
	Q_EQ   = iota // quartets equal
	Q_NEQ         // quartets not equal
	Q_DIFF        // quartets on different taxa set
)

/*
	create new quartet

args:

	*tree.Tree qTree (quartet tree)
	*tree.Tree tre   (reference tree with taxa set)
*/
func NewQuartet(qTree, tre *tree.Tree) (*Quartet, error) {
	qTaxa := qTree.AllTipNames()
	if len(qTaxa) != 4 {
		return nil, fmt.Errorf("Bad quartet, tree has %d != 4 leaves", len(qTaxa))
	}
	qTree.UnRoot()
	if qTree.Root().Nneigh() != 3 {
		return nil, errors.New("Bad quartet, probably does not contain bipartition")
	}
	taxaIDs := make([]int, 4)
	leaves := qTree.Tips()
	idToBool := make(map[int]bool) // true is one side of the bipartition, false is the other
	for i, l := range leaves {
		ti, err := tre.TipIndex(l.Name())
		tipIndexPanic("convertQuartet", err)
		taxaIDs[i] = ti
		r, err := l.Parent()
		if err != nil && err.Error() == "The node has more than one parent" { // we ignore the error produced when cur = root
			panic(fmt.Errorf("convertQuartet: %w", err))
		}
		idToBool[ti] = r == qTree.Root()
	}
	slices.Sort(taxaIDs) // sort ids so quartet topologies are equal if they are the same
	count := 0           // only to varify we have exactly two taxa on each side
	var power uint8 = 1
	var topo uint8 = TOPONIL
	for _, id := range taxaIDs {
		if idToBool[id] {
			count++
			topo |= power
		}
		power *= 2
	}
	if count != 2 {
		panic("quartet didn't define bipartition properly, probably due to a bug")
	}
	if topo != TOPO1 && topo != TOPO2 && topo != TOPO3 {
		topo ^= 0b1111 // inverse so that topology is consistent with enumerated constants
	}

	return &Quartet{taxa: taxaIDs, topology: topo}, nil
}

func tipIndexPanic(context string, err error) {
	if err != nil {
		panic(fmt.Errorf("%s TipIndex panic: %w", context, err))
	}
}

/* Not efficent, do no use except for testing !!! */
func (q *Quartet) String(tre *tree.Tree) string {
	names := make(map[int]string)
	for _, u := range tre.Tips() {
		ti, err := tre.TipIndex(u.Name())
		tipIndexPanic("String", err)
		names[ti] = u.Name()
	}
	qString := "|"
	for i := 0; i < 4; i++ {
		if (q.topology>>i)%2 == 0 {
			qString += names[q.taxa[i]]
		} else {
			qString = names[q.taxa[i]] + qString
		}
	}
	return qString
}

/*
	 compares two quartets; there are three possible results:
		- Q_DIFF (they contain different taxa)
		- Q_NEQ  (they have a different topology)
		- Q_EQ   (they have the same topology)
*/
func (q1 *Quartet) Compare(q2 *Quartet) int {
	for i := 0; i < 4; i++ {
		if q1.taxa[i] != q2.taxa[i] {
			return Q_DIFF
		}
	}
	if q1.topology^q2.topology != 0b1111 && q1.topology^q2.topology != 0b0000 {
		return Q_NEQ
	} else {
		return Q_EQ
	}
}
