package score

import (
	"golang.org/x/sync/errgroup"

	gr "github.com/jsdoublel/camus/internal/graphs"
)

func CalcuateEdgePenalties(td *gr.TreeData, nprocs int) ([][]uint, error) {
	n := len(td.Nodes())
	edgePenalties := make([][]uint, n)
	var g errgroup.Group
	g.SetLimit(nprocs)
	for u := range n {
		g.Go(func() error {
			edgePenalties[u] = make([]uint, n)
			for w := range n {
				if shouldCalcEdge(u, w, td) {
					edgePenalties[u][w] = calculatePenalty(u, w, td)
				}
			}
			return nil
		})
	}
	return edgePenalties, g.Wait()
}

// Calculates the number of quartets that *could* be added by the addition of
// this edge (from u to w)..
func calculatePenalty(u, w int, td *gr.TreeData) uint {
	subsets := getNumTaxaUnderNodes(u, w, td)
	if len(subsets) < 4 {
		panic("less than four subsets; this should not happen")
	}
	coe := [...]uint{1, 0, 0, 0}
	for i := 1; i < len(subsets); i++ {
		for j := 3; j > 0; j-- {
			coe[j] += subsets[i] * coe[j-1]
		}
	}
	return subsets[0] * coe[3]
}

// Gets the number of nodes in each subtree connected to the unrooted cycle
// formed by the edge. Nodes under w will be at index 0.
func getNumTaxaUnderNodes(u, w int, td *gr.TreeData) (subsets []uint) {
	v := td.LCA(u, w)
	if u == v {
		parent, err := td.IdToNodes[v].Parent()
		if err != nil {
			panic(err)
		}
		subsets = append(collectNodesUpPath(w, parent.Id(), td), uint(td.NLeaves)-td.NumLeavesBelow[v])
		return
	}
	subsets = append(collectNodesUpPath(w, v, td), collectNodesUpPath(u, v, td)...)
	subsets = append(subsets, uint(td.NLeaves)-td.NumLeavesBelow[v])
	return
}

func collectNodesUpPath(start, end int, td *gr.TreeData) (subsets []uint) {
	subsets = []uint{td.NumLeavesBelow[start]}
	cur := start
	parentID := -1
	for parentID != end {
		parent, err := td.IdToNodes[cur].Parent()
		if err != nil {
			if err.Error() == "The node has more than one parent" {
				panic(err)
			}
			panic("end is below start while adding up nodes on path")
		}
		parentID = parent.Id()
		c := td.Children[parentID]
		if c[0].Id() == cur {
			subsets = append(subsets, td.NumLeavesBelow[c[1].Id()])
		} else {
			subsets = append(subsets, td.NumLeavesBelow[c[0].Id()])
		}
		cur = parentID
	}
	return
}
