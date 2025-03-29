package qrt

import (
	"errors"
	"fmt"

	"github.com/evolbioinfo/gotree/tree"
)

type Quartet struct {
	Taxa     [NTaxa]int // should be in sorted order
	Topology uint8      // represents one of three possible quartet topologies
}

const (
	NTaxa = 4
	NTopo = 3

	// Possible results for quartet comparison (currently unused)
	Qeq   = iota // quartets equal
	Qneq         // quartets not equal
	Qdiff        // quartets on different taxa set

	Qtopo1 = uint8(0b1100)
	Qtopo2 = uint8(0b1010)
	Qtopo3 = uint8(0b0110)
)

var (
	ErrTipNameMismatch = errors.New("tip name mismatch! maybe the gene tree and constraint tree labels don't match?")
	ErrInvalidQuartet  = errors.New("invalid newick for quartet")
)

// Generates quartet from four leaf newick tree (only used for testing)
func NewQuartet(qTree, tre *tree.Tree) (*Quartet, error) {
	qTaxa := qTree.AllTipNames()
	if len(qTaxa) != NTaxa {
		return nil, fmt.Errorf("%w, tree has %d != 4 leaves", ErrInvalidQuartet, len(qTaxa))
	}
	qTree.UnRoot()
	if qTree.Root().Nneigh() != 3 {
		return nil, fmt.Errorf("%w, probably does not contain bipartition", ErrInvalidQuartet)
	}
	taxaIDs := [NTaxa]int{}
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
	topo := setTopology(&taxaIDs)
	return &Quartet{Taxa: taxaIDs, Topology: topo}, nil
}

// Generate unit8 representing quartet topology
func setTopology(taxaIDs *[NTaxa]int) uint8 {
	if len(taxaIDs) != NTaxa {
		panic("taxaIDs len != 4 in setTopology")
	}
	topo := sortTaxa(taxaIDs) // sort ids so quartet topologies are equal if they are the same
	if topo%2 != 0 {          // normalize quartet (i.e., so that there are three topologies instead of six)
		topo ^= 0b1111
	}
	if topo != Qtopo1 && topo != Qtopo2 && topo != Qtopo3 {
		panic(fmt.Sprintf("quartet didn't define bipartition properly, probably due to a bug: %b", topo))
	}
	return topo
}

// Short 4 long int array (no build in array sort in go)
// returns the topology as uint8
func sortTaxa(arr *[NTaxa]int) uint8 {
	topo := uint8(0b0011)
	for i := 0; i < 3; i++ {
		for j := i + 1; j < 4; j++ {
			if arr[i] > arr[j] {
				if topo>>i&1 != topo>>j&1 {
					topo ^= (1 << i) | (1 << j)
				}
				arr[i], arr[j] = arr[j], arr[i]
			}
		}
	}
	return topo
}

// Returns hashmap containing quartets from tree
func QuartetsFromTree(tre, constTree *tree.Tree) (map[Quartet]uint, error) {
	tre.UnRoot()                           // some quartets are missed if tree is rooted
	treeQuartets := make(map[Quartet]uint) // get quartets from tree
	taxaIDsMap, err := mapIDsFromConstTree(tre, constTree)
	if err != nil {
		return nil, err
	}
	tre.Quartets(false, func(q *tree.Quartet) {
		treeQuartets[*quartetFromTreeQ(q, taxaIDsMap)] = 1
	})
	return treeQuartets, nil
}

// Create quartet from gotree *tree.Quartet
func quartetFromTreeQ(tq *tree.Quartet, constMap []int) *Quartet {
	taxaIDs := [...]int{constMap[int(tq.T1)], constMap[int(tq.T2)], constMap[int(tq.T3)], constMap[int(tq.T4)]}
	return &Quartet{Taxa: taxaIDs, Topology: setTopology(&taxaIDs)}
}

func mapIDsFromConstTree(gtre, tre *tree.Tree) ([]int, error) {
	nLeavesGtree, err := gtre.NbTips()
	if err != nil {
		panic(fmt.Sprintf("gene tree %s", err))
	}
	idMap := make([]int, nLeavesGtree)
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

// Not efficent, do no use except for testing !!!
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
	for i := range NTaxa {
		if (q.Topology>>i)%2 == 0 {
			qString += names[q.Taxa[i]]
		} else {
			qString = names[q.Taxa[i]] + qString
		}
	}
	return qString
}

func QSetToString(qSet map[Quartet]uint, tre *tree.Tree) string {
	str := "{"
	for q, c := range qSet {
		str += fmt.Sprintf("%s:%d, ", q.String(tre), c)
	}
	return str[:len(str)-2] + "}"
}

// Compares two quartets (currently unused).
// There are three possible results:
//   - Qdiff (they contain different taxa)
//   - Qneq  (they have a different topology)
//   - Qeq   (they have the same topology)
func (q1 *Quartet) Compare(q2 *Quartet) int {
	for i := range NTaxa {
		if q1.Taxa[i] != q2.Taxa[i] {
			return Qdiff
		}
	}
	if q1.Topology^q2.Topology != 0b1111 && q1.Topology^q2.Topology != 0b0000 {
		return Qneq
	} else {
		return Qeq
	}
}
