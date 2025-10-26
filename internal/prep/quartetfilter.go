package prep

import (
	"cmp"
	"fmt"
	"slices"
	"strconv"

	gr "github.com/jsdoublel/camus/internal/graphs"
)

// Options for quartet filter mode
type QuartetFilterOptions struct {
	mode      QMode     // mode (value between 0 and 1)
	threshold Threshold // threshold for filtering [0, 1]
}

func SetQuartetFilterOptions(mode int, threshold float64) (QuartetFilterOptions, error) {
	var m QMode
	if err := m.Set(mode); err != nil {
		return QuartetFilterOptions{}, err
	}
	var t Threshold
	if err := t.Set(threshold); err != nil {
		return QuartetFilterOptions{}, err
	}
	return QuartetFilterOptions{mode: m, threshold: t}, nil
}

func (opts QuartetFilterOptions) QuartetFilterOff() bool {
	return opts.mode == 0
}

type QMode int

const (
	NonRestrictive QMode = iota + 1
	Restrictive
)

func (mode *QMode) Set(n int) error {
	if n < 0 || n > 2 {
		return fmt.Errorf("quartet mode %d is %w", n, ErrTypeOutRange)
	}
	*mode = QMode(n)
	return nil
}

func (mode QMode) String() string {
	return strconv.Itoa(int(mode))
}

type Threshold float64

func (thresh *Threshold) Set(n float64) error {
	if n < 0 || n > 1 {
		return fmt.Errorf("threshold %f is %w", n, ErrTypeOutRange)
	}
	*thresh = Threshold(n)
	return nil
}

func (thresh Threshold) String() string {
	return strconv.FormatFloat(float64(thresh), 'f', -1, 64)
}

func (thresh Threshold) Keep(counts []uint32) bool {
	if len(counts) != 3 {
		panic("there should be three counts, one for each quartet topology")
	}
	slices.Sort(counts)
	sum := counts[0] + counts[1]
	return uint32(float64(thresh)*float64(sum)) < counts[1]-counts[0]
}

func filterQuartets(qCounts map[gr.Quartet]uint32, opts QuartetFilterOptions) {
	for q := range qCounts {
		quartets := q.AllQuartets()
		counts := []uint32{qCounts[quartets[0]], qCounts[quartets[1]], qCounts[quartets[2]]}
		slices.SortFunc(quartets, func(q1, q2 gr.Quartet) int {
			return cmp.Compare(qCounts[q1], qCounts[q2])
		})
		if !opts.threshold.Keep(counts) {
			delete(qCounts, quartets[0])
			delete(qCounts, quartets[1])
			continue
		}
		switch opts.mode {
		case NonRestrictive:
		case Restrictive:
			delete(qCounts, quartets[0])
		default:
			panic("invalid quartet mode case")
		}
	}
}
