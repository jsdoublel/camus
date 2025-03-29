package qrt

import (
	"errors"
	"fmt"
)

var (
	mapper     map[uint8]uint8 = map[uint8]uint8{Qtopo1: 0, Qtopo2: 1, Qtopo3: 2}
	maxSafeInt                 = uint(1 << 53)

	ErrCannotConvert = errors.New("cannot convert")
)

// wrapper for recording quartet count stats
type QuartetStats struct {
	qstats map[[NTaxa]int]*[NTopo]uint // map mapping hash of quartet taxa to list giving count of the three topologies
	totals map[[NTaxa]int]uint
}

// Makes empty QuartetStats struct
func MakeQStats() *QuartetStats {
	return &QuartetStats{qstats: make(map[[NTaxa]int]*[NTopo]uint), totals: make(map[[NTaxa]int]uint)}
}

// Add to count of particular quartet
func (qs *QuartetStats) Count(q *Quartet, count uint) {
	if qs.qstats[q.Taxa] == nil {
		qs.qstats[q.Taxa] = new([NTopo]uint)
	}
	qs.qstats[q.Taxa][mapper[q.Topology]] += count
	qs.totals[q.Taxa] += count
}

// Get count of particular quartet
func (qs *QuartetStats) GetCount(q *Quartet) uint {
	return qs.qstats[q.Taxa][mapper[q.Topology]]
}

func (qs *QuartetStats) Percent(q *Quartet) float64 {
	total := uint(0)
	for k := range mapper {
		total += qs.qstats[q.Taxa][mapper[k]]
	}
	fTotal, err := uintToFloat64Safe(total)
	if err != nil {
		panic(err)
	}
	fNumerator, err := uintToFloat64Safe(qs.qstats[q.Taxa][mapper[q.Topology]])
	if err != nil {
		panic(err)
	}
	result := fNumerator / fTotal
	if result < 0 || result > 1 {
		panic(fmt.Sprintf("percent value %f out of range", result))
	}
	return result
}

func uintToFloat64Safe(n uint) (float64, error) {
	if n >= maxSafeInt {
		return 0, fmt.Errorf("%w %d to float64 exactly", ErrCannotConvert, n)
	}
	return float64(n), nil
}
