package score

import (
	"context"

	"github.com/evolbioinfo/gotree/tree"
	"golang.org/x/sync/errgroup"

	gr "github.com/jsdoublel/camus/internal/graphs"
)

const Max16Bit = ^uint16(0)

// Calculate scores for all edges
func CalculateEdgeScores(td *gr.TreeData, nprocs int) ([][]uint, error) {
	n := len(td.Nodes())
	edgeScores := make([][]uint, n)
	g, _ := errgroup.WithContext(context.Background())
	g.SetLimit(nprocs)
	for u := range n {
		g.Go(func() error {
			edgeScores[u] = make([]uint, n)
			for w := range n {
				if shouldCalcEdge(u, w, td) {
					edgeScores[u][w] = getEdgeScore(u, w, td)
				}
			}
			return nil
		})
	}
	return edgeScores, g.Wait()
}

func shouldCalcEdge(u, w int, td *gr.TreeData) bool {
	return !td.Under(w, u) && CycleLength(u, w, td) > 3
}

func CycleLength(u, w int, td *gr.TreeData) int {
	v := td.LCA(u, w)
	length := (td.Depths[u] - td.Depths[v]) + (td.Depths[w] - td.Depths[v]) + 1
	if v == u { // we have to account for the edge above v that our new edge is anchored to
		length += 1
	}
	return length
}

func getEdgeScore(u, w int, td *gr.TreeData) uint {
	v := td.LCA(u, w)
	uNode, wNode, vNode := td.IdToNodes[u], td.IdToNodes[w], td.IdToNodes[v]
	edgeScore := uint(0)
	wSub := getWSubtree(u, w, v, td)
	for _, q := range td.Quartets(v) {
		if QuartetScore(q, uNode, wNode, vNode, wSub, td) == gr.Qeq {
			edgeScore += td.NumQuartet(q)
		}
	}
	return edgeScore
}

func getWSubtree(u, w, v int, td *gr.TreeData) *tree.Node {
	switch {
	case u == v:
		return td.IdToNodes[v]
	case td.Under(td.Children[v][0].Id(), w) || w == td.Children[v][0].Id():
		return td.IdToNodes[td.Children[v][0].Id()]
	default:
		return td.IdToNodes[td.Children[v][1].Id()]
	}
}

func QuartetScore(q gr.Quartet, u, w, v, wSub *tree.Node, td *gr.TreeData) int {
	bottom, bi, unique := uniqueTaxaBelowNodeFromQ(w, q, td)
	if !unique || bottom == Max16Bit {
		return gr.Qdiff
	}
	cycleNodes := make(map[int]bool)
	taxaToLCA := make(map[uint16]int) // tip index -> lca
	for _, t := range q.Taxa() {
		tID := td.NodeID(t)
		var lca int
		switch {
		case !td.InLeafset(uint16(v.Id()), t):
			lca = 0
		case td.InLeafset(uint16(wSub.Id()), t) || td.InLeafset(uint16(u.Id()), uint16(bottom)):
			lca = td.LCA(w.Id(), tID)
		default:
			lca = td.LCA(u.Id(), tID)
		}
		cycleNodes[lca] = true
		taxaToLCA[t] = lca
	}
	if len(cycleNodes) != 4 {
		return gr.Qdiff
	}
	neighbor := neighborTaxaQ(q, bi)
	lcaDepths := make(map[int]int) // node ID -> depth
	for k, v := range cycleNodes {
		if v {
			lcaDepths[k] = td.Depths[k]
		}
	}
	nLeaves := td.NLeaves
	minW, maxU := nLeaves, -1
	var bestTaxa uint16
	taxaInU := false
	for _, t := range q.Taxa() {
		d := lcaDepths[taxaToLCA[t]]
		if !taxaInU && (td.InLeafset(uint16(wSub.Id()), t) && d < minW) {
			minW = d
			bestTaxa = t
		} else if !td.InLeafset(uint16(wSub.Id()), t) && d > maxU {
			taxaInU = true
			maxU = d
			bestTaxa = t
		}
	}
	if bestTaxa == neighbor {
		return gr.Qeq
	} else {
		return gr.Qneq
	}
}

// Returns -1 for both id and index if no taxa is found, true if taxa is unique (or there isn't a taxa)
func uniqueTaxaBelowNodeFromQ(n *tree.Node, q gr.Quartet, td *gr.TreeData) (uint16, int, bool) {
	taxaID := Max16Bit
	taxaIndex := -1
	for i, t := range q.Taxa() {
		if td.InLeafset(uint16(n.Id()), t) && taxaID == Max16Bit {
			taxaID, taxaIndex = t, i
		} else if td.InLeafset(uint16(n.Id()), t) {
			return taxaID, taxaIndex, false
		}
	}
	return taxaID, taxaIndex, true
}

// Return neighbor of taxa at index i in quartet
func neighborTaxaQ(q gr.Quartet, i int) uint16 {
	b := (q.Topology() >> i) % 2
	for j := range 4 {
		if j != i && (q.Topology()>>j)%2 == b {
			return q.Taxon(j)
		}
	}
	panic("invalid quartet or bad i")
}
