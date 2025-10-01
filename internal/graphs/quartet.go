package graphs

import (
	"errors"
	"fmt"
	"iter"
	"log"

	"github.com/evolbioinfo/gotree/tree"
)

type Quartet uint64

const (
	NilQuartet = 0
	NTaxa      = 4

	taxaShift = 15
	topoShift = 60

	taxaMask = (1 << taxaShift) - 1 // 0x7FFF
	topoMask = (1 << topoShift) - 1 // 0xF

	Qtopo1 = uint8(0b1100) // three quartet topologies
	Qtopo2 = uint8(0b1010)
	Qtopo3 = uint8(0b0110)

	// Possible results for quartet comparison
	Qneq  = iota // quartets equal
	Qeq          // quartets not equal
	Qdiff        // quartets on different taxa set
)

var (
	ErrTipNameMismatch = errors.New("tip name mismatch! maybe the gene tree and constraint tree labels don't match?")
	ErrInvalidQuartet  = errors.New("invalid newick for quartet")
)

// Generates quartet from four leaf newick tree (only used for testing)
func NewQuartet(qTree, tre *tree.Tree) (Quartet, error) {
	qTaxa := qTree.AllTipNames()
	if len(qTaxa) != 4 {
		return 0, fmt.Errorf("%w, tree has %d != 4 leaves", ErrInvalidQuartet, len(qTaxa))
	}
	qTree.UnRoot()
	if qTree.Root().Nneigh() != 3 {
		return 0, fmt.Errorf("%w, probably does not contain bipartition", ErrInvalidQuartet)
	}
	taxaIDs := [4]int16{}
	leaves := qTree.Tips()
	idToBool := make(map[int]bool) // true is one side of the bipartition, false is the other
	for i, l := range leaves {
		ti, err := tre.TipIndex(l.Name())
		if err != nil {
			panic(err)
		}
		taxaIDs[i] = int16(ti)
		r, err := l.Parent()
		if err != nil && err.Error() == "The node has more than one parent" { // we ignore the error produced when cur = root
			panic(fmt.Errorf("convertQuartet: %w", err))
		}
		idToBool[ti] = r == qTree.Root()
	}
	topo := setTopology(&taxaIDs)
	return makeQuartet(taxaIDs, topo), nil
}

func makeQuartet(taxa [4]int16, topology uint8) Quartet {
	var q uint64
	for i, t := range taxa {
		q |= uint64(t) << (taxaShift * i) // we assume positive taxa ids
	}
	q |= uint64(topology) << topoShift
	return Quartet(q)
}

// Generate unit8 representing quartet topology
func setTopology(taxaIDs *[4]int16) uint8 {
	if len(taxaIDs) != 4 {
		panic("taxaIDs len != 4 in setTopology")
	}
	topo := sortTaxa(taxaIDs) // sort ids so quartet topologies are equal if they are the same
	if topo%2 != 0 {          // normalize quartet (i.e., so that there are three topologies instead of six)
		topo ^= 0b1111
	}
	if topo != 0b1100 && topo != 0b0110 && topo != 0b1010 {
		panic(fmt.Sprintf("quartet didn't define bipartition properly, probably due to a bug: %b", topo))
	}
	return topo
}

// Short 4 long int array (no build in array sort in go)
// returns the topology as uint8
func sortTaxa(arr *[4]int16) uint8 {
	topo := uint8(0b0011)
	for i := 0; i < 3; i++ {
		for j := i + 1; j < 4; j++ {
			if arr[i] > arr[j] {
				bi := uint8(topo >> i & 1)
				bj := uint8(topo >> j & 1)
				if bi != bj {
					m := uint8((1 << i) | (1 << j))
					topo ^= m
				}
				arr[i], arr[j] = arr[j], arr[i]
			}
		}
	}
	return topo
}

func missmatchTaxaSets(tre1, tre2 *tree.Tree) (bool, error) {
	n1, err := tre1.NbTips()
	if err != nil {
		return false, err
	}
	n2, err := tre2.NbTips()
	if err != nil {
		return false, err
	}
	return n1 != n2, nil
}

// Returns hashmap containing quartets from tree
func QuartetsFromTree(tre, constTree *tree.Tree) (map[Quartet]uint32, error) {
	if b, err := missmatchTaxaSets(tre, constTree); err != nil {
		return nil, err
	} else if b {
		log.Println("WARNING: missing taxa detected in one or more gene trees; this may cause issues with some scoring metrics")
	}
	tre.UnRoot() // some quartets are missed if tree is rooted
	treeQuartets := make(map[Quartet]uint32)
	taxaIDsMap, err := MapIDsFromConstTree(tre, constTree)
	if err != nil {
		return nil, err
	}
	tre.Quartets(false, func(q *tree.Quartet) {
		treeQuartets[QuartetFromTreeQ(q, taxaIDsMap)] = 1
	})
	return treeQuartets, nil
}

// Create quartet from gotree *tree.Quartet
func QuartetFromTreeQ(tq *tree.Quartet, constMap []int16) Quartet {
	taxaIDs := [...]int16{constMap[tq.T1], constMap[tq.T2], constMap[tq.T3], constMap[tq.T4]}
	return makeQuartet(taxaIDs, setTopology(&taxaIDs))
}

func MapIDsFromConstTree(gtre, tre *tree.Tree) ([]int16, error) {
	nLeavesGtree, err := gtre.NbTips()
	if err != nil {
		panic(fmt.Sprintf("gene tree %s", err))
	}
	idMap := make([]int16, nLeavesGtree)
	for _, name := range gtre.AllTipNames() {
		constTreeID, err := tre.TipIndex(name)
		if err != nil {
			return nil, fmt.Errorf("%w, %s", ErrTipNameMismatch, err.Error())
		}
		gTreeID, err := gtre.TipIndex(name)
		idMap[gTreeID] = int16(constTreeID)
		if err != nil {
			return nil, fmt.Errorf("%w, %s", ErrTipNameMismatch, err.Error())
		}
	}
	return idMap, nil
}

func (q Quartet) Topology() uint8 {
	return uint8((q >> topoShift) & topoMask)
}

func (q Quartet) Taxon(i int) uint16 {
	return uint16((q >> (taxaShift * i)) & taxaMask)
}

func (q Quartet) Taxa() iter.Seq2[int, uint16] {
	return func(yield func(int, uint16) bool) {
		for i := range 4 {
			if !yield(i, q.Taxon(i)) {
				return
			}
		}
	}
}

func (q Quartet) AllQuartets() []Quartet {
	// Use bit operations: keep taxa bits, replace topology bits.
	base := uint64(q) & ^(uint64(0xF) << topoShift)
	return []Quartet{
		Quartet(base | (uint64(Qtopo1) << topoShift)),
		Quartet(base | (uint64(Qtopo2) << topoShift)),
		Quartet(base | (uint64(Qtopo3) << topoShift)),
	}
}

// Not efficient, do no use except for testing !!!
func (q *Quartet) String(tre *tree.Tree) string {
	names := make(map[uint16]string)
	for _, u := range tre.Tips() {
		ti, err := tre.TipIndex(u.Name())
		if err != nil {
			panic(err)
		}
		names[uint16(ti)] = u.Name()
	}
	qString := "|"
	for i := range 4 {
		if (q.Topology()>>i)%2 == 0 {
			qString += names[q.Taxon(i)]
		} else {
			qString = names[q.Taxon(i)] + qString
		}
	}
	return qString
}

func QSetToString(qSet map[Quartet]uint32, tre *tree.Tree) string {
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
func (q1 *Quartet) Compare(q2 Quartet) int {
	for i := range 4 {
		if q1.Taxon(i) != q2.Taxon(i) {
			return Qdiff
		}
	}
	if q1.Topology()^q2.Topology() != 0b1111 && q1.Topology()^q2.Topology() != 0b0000 {
		return Qneq
	} else {
		return Qeq
	}
}
