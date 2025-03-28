package prep

import (
	"errors"
	"fmt"

	"github.com/evolbioinfo/gotree/tree"
)

type Quartet struct {
	Taxa     [4]int // should be in sorted order
	Topology uint8  // represents one of three possible quartet topologies
}

// Possible results for quartet comparison (currently unused)
const (
	Q_EQ   = iota // quartets equal
	Q_NEQ         // quartets not equal
	Q_DIFF        // quartets on different taxa set
)

var (
	ErrTipNameMismatch = errors.New("tip name mismatch! maybe the gene tree and constraint tree labels don't match?")
	ErrInvalidQuartet  = errors.New("invalid newick for quartet")
)

// Generates quartet from four leaf newick tree (only used for testing)
func NewQuartet(qTree, tre *tree.Tree) (*Quartet, error) {
	qTaxa := qTree.AllTipNames()
	if len(qTaxa) != 4 {
		return nil, fmt.Errorf("%w, tree has %d != 4 leaves", ErrInvalidQuartet, len(qTaxa))
	}
	qTree.UnRoot()
	if qTree.Root().Nneigh() != 3 {
		return nil, fmt.Errorf("%w, probably does not contain bipartition", ErrInvalidQuartet)
	}
	taxaIDs := [4]int{}
	leaves := qTree.Tips()
	idToBool := make(map[int]bool) // true is one side of the bipartition, false is the other
	for i, l := range leaves {
		ti, err := tre.TipIndex(l.Name())
		if err != nil {
			panic(err)
		}
		taxaIDs[i] = ti
		r, err := l.Parent()
		if err != nil && err.Error() == "The node has more than one parent" { // we ignore the error produced when cur = root
			panic(fmt.Errorf("convertQuartet: %w", err))
		}
		idToBool[ti] = r == qTree.Root()
	}
	topo := setTopology(&taxaIDs, idToBool)
	return &Quartet{Taxa: taxaIDs, Topology: topo}, nil
}

// Generate unit8 representing quartet topology
func setTopology(taxaIDs *[4]int, idToBool map[int]bool) uint8 {
	if len(taxaIDs) != 4 {
		panic("taxaIDs len != 4 in setTopology")
	}
	sortTaxa(taxaIDs) // sort ids so quartet topologies are equal if they are the same
	count := 0        // only to varify we have exactly two taxa on each side
	var power uint8 = 1
	var topo uint8 = 0b0000
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
	if topo%2 != 0 { // normalize quartet (i.e., so that there are three topologies instead of six)
		topo ^= 0b1111
	}
	return topo
}

// Short 4 long int array (no build in array sort in go)
func sortTaxa(arr *[4]int) {
	for i := 0; i < 3; i++ {
		for j := i + 1; j < 4; j++ {
			if arr[i] > arr[j] {
				arr[i], arr[j] = arr[j], arr[i]
			}
		}
	}
}

// Returns hashmap containing quartets from tree
func QuartetsFromTree(tre, constTree *tree.Tree) (map[Quartet]uint, error) {
	tre.UnRoot()                           // some quartets are missed if tree is rooted
	treeQuartets := make(map[Quartet]uint) // get quartets from tree
	taxaIDsMap, err := mapIDsFromConstTree(tre, constTree)
	if err != nil {
		return nil, err
	}
	// n, err := tre.NbTips()
	// if err != nil {
	// 	panic(err)
	// }
	// leafsets := make([][]int, n)
	// tre.PostOrder(func(cur, prev *tree.Node, e *tree.Edge) (keep bool) {
	// 	leafsets[cur.Id()] = make([]int, 0)
	// 	if cur.Tip() {
	// 		leafsets[cur.Id()] = append(leafsets[cur.Id()], cur.Id())
	// 	} else if cur != tre.Root() {
	// 		leafsets[cur.Id()] = append(leafsets[cur.Id()], )
	// 	}
	// 	return true
	// })
	tre.Quartets(false, func(q *tree.Quartet) {
		treeQuartets[*quartetFromTreeQ(q, taxaIDsMap)] = 1
	})
	return treeQuartets, nil
}

// Create quartet from gotree *tree.Quartet
func quartetFromTreeQ(tq *tree.Quartet, constMap map[int]int) *Quartet {
	taxaIDs := [...]int{constMap[int(tq.T1)], constMap[int(tq.T2)], constMap[int(tq.T3)], constMap[int(tq.T4)]}
	idToBool := make(map[int]bool)
	idToBool[taxaIDs[0]] = true
	idToBool[taxaIDs[1]] = true
	return &Quartet{Taxa: taxaIDs, Topology: setTopology(&taxaIDs, idToBool)}
}

func mapIDsFromConstTree(gtre, tre *tree.Tree) (map[int]int, error) {
	idMap := make(map[int]int)
	for _, name := range gtre.AllTipNames() {
		constTreeID, err := tre.TipIndex(name)
		if err != nil {
			return nil, fmt.Errorf("%w, %s", ErrTipNameMismatch, err.Error())
		}
		gTreeID, err := gtre.TipIndex(name)
		idMap[gTreeID] = constTreeID
		if err != nil {
			return nil, fmt.Errorf("%w, %s", ErrTipNameMismatch, err.Error())
		}
	}
	return idMap, nil
}

/* Not efficent, do no use except for testing !!! */
func (q *Quartet) String(tre *tree.Tree) string {
	names := make(map[int]string)
	for _, u := range tre.Tips() {
		ti, err := tre.TipIndex(u.Name())
		if err != nil {
			panic(err)
		}
		names[ti] = u.Name()
	}
	qString := "|"
	for i := range 4 {
		if (q.Topology>>i)%2 == 0 {
			qString += names[q.Taxa[i]]
		} else {
			qString = names[q.Taxa[i]] + qString
		}
	}
	return qString
}

// Compares two quartets (currently unused).
// There are three possible results:
//   - Q_DIFF (they contain different taxa)
//   - Q_NEQ  (they have a different topology)
//   - Q_EQ   (they have the same topology)
func (q1 *Quartet) Compare(q2 *Quartet) int {
	for i := range 4 {
		if q1.Taxa[i] != q2.Taxa[i] {
			return Q_DIFF
		}
	}
	if q1.Topology^q2.Topology != 0b1111 && q1.Topology^q2.Topology != 0b0000 {
		return Q_NEQ
	} else {
		return Q_EQ
	}
}
