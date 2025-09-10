package score

import (
	"context"

	"golang.org/x/sync/errgroup"

	gr "github.com/jsdoublel/camus/internal/graphs"
)

func CalcuateEdgePenalties(td *gr.TreeData, nprocs int) ([][]uint, error) {
	n := len(td.Nodes())
	edgePenalties := make([][]uint, n)
	g, _ := errgroup.WithContext(context.Background())
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

func calculatePenalty(u, w int, td *gr.TreeData) uint {
	return 0
}

func getNumTaxaUnderNodes(u, w int, td *gr.TreeData) []uint {
	v := td.LCA(u, w)
	result := append(addNodesUpPath(u, v, td), addNodesUpPath(w, v, td)...)
	return append(result, td.NumLeavesBelow[v])
}

func addNodesUpPath(start, end int, td *gr.TreeData) []uint {
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
