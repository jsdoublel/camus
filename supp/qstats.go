package supp

import "github.com/jsdoublel/camus/prep"

var mapper map[uint8]uint8 = map[uint8]uint8{0b1010: 0, 0b1100: 1, 0b0110: 2}

// wrapper for recording quartet count stats
type QuartetStats struct {
	qstats map[[4]int][3]uint // map mapping hash of quartet taxa to list giving count of the three topologies
}

// Add to count of particular quartet
func (qs *QuartetStats) Count(q *prep.Quartet, count uint) {
	counter := qs.qstats[q.Taxa]
	counter[mapper[q.Topology]] += count
}

// Get count of particular quartet
func (qs *QuartetStats) GetCount(q *prep.Quartet) uint {
	return qs.qstats[q.Taxa][mapper[q.Topology]]
}
