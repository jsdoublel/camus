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
func getNumTaxaUnderNodes(u, w int, td *gr.TreeData) []uint {
	v := td.LCA(u, w)
	result := append(collectNodesUpPath(w, v, td), collectNodesUpPath(u, v, td)...)
	return append(result, td.NumLeavesBelow[v])
}

func collectNodesUpPath(start, end int, td *gr.TreeData) []uint {
	result := make([]uint, 1)
	result[0] = td.NumLeavesBelow[start]
	cur := start
	for cur != end {
		p, err := td.IdToNodes[cur].Parent()
		if err.Error() == "The node has more than one parent" {
			panic(err)
		} else if err != nil {
			panic("end is below start while adding up nodes on path")
		}
		c := td.Children[p.Id()]
		if c[0].Id() == cur {
			result = append(result, td.NumLeavesBelow[c[1].Id()])
		} else {
			result = append(result, td.NumLeavesBelow[c[0].Id()])
		}
	}
	return result
}
