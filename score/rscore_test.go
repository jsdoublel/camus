package score

import (
	"math"
	"strings"
	"testing"

	"github.com/evolbioinfo/gotree/io/newick"
	"github.com/evolbioinfo/gotree/tree"

	"github.com/jsdoublel/camus/prep"
)

func TestCalculateRecticulationScore(t *testing.T) {
	testCases := []struct {
		name     string
		network  string
		gtrees   []string
		expected []*map[string]float64
	}{
		{
			name:    "basic",
			network: "(((9,0),(7,(6,(#H0,8h0u)))),((#H2,(12,((3,(14h2w)#H2),10))h2u),((((5,(#H1,13h1u)),((2h1w)#H1,11))h0w)#H0,(1,4))));",
			gtrees: []string{
				"((0,9),(5,7));",
				"((5,7),(9,6));",
				"((9,7),(5,6));",
				"((5,9),(7,6));",
			},
			expected: []*map[string]float64{
				{"#H0": math.NaN(), "#H1": math.NaN(), "#H2": math.NaN()},
				{"#H0": float64(0), "#H1": math.NaN(), "#H2": math.NaN()},
				{"#H0": float64(1), "#H1": math.NaN(), "#H2": math.NaN()},
				{"#H0": float64(0), "#H1": math.NaN(), "#H2": math.NaN()},
			},
		},
	}
	for _, test := range testCases {
		t.Run(test.name, func(t *testing.T) {
			tre, err := newick.NewParser(strings.NewReader(test.network)).Parse()
			if err != nil {
				t.Fatalf("invalid newick in file %s", err)
			}
			ntw, err := prep.ConvertToNetwork(tre)
			if err != nil {
				t.Fatalf("test case failed with unexpected error %s", err)
			}
			gtrees := make([]*tree.Tree, len(test.gtrees))
			for i, gt := range test.gtrees {
				tmp, err := newick.NewParser(strings.NewReader(gt)).Parse()
				if err != nil {
					t.Fatal("invalid newick tree; test is written wrong")
				}
				gtrees[i] = tmp
			}
			result, err := CalculateReticulationScore(ntw, gtrees)
			if err != nil {
				t.Errorf("test case failed with unexpected error %s", err)
			}
			if !compareScoreMaps(result, test.expected) {
				t.Error("result != expected", result, test.expected)
			}

		})
	}
}

// compares the two maps (specifically allows NaN == NaN to be true)
func compareScoreMaps(m1, m2 []*map[string]float64) bool {
	if len(m1) != len(m2) {
		return false
	}
	for i := range m1 {
		if len(*m1[i]) != len(*m2[i]) {
			return false
		}
		for key := range *m1[i] {
			if (*m1[i])[key] != (*m2[i])[key] && (!math.IsNaN((*m1[i])[key]) || !math.IsNaN((*m2[i])[key])) {
				return false
			}
		}
	}
	return true
}

func BenchmarkCalculateRecticulationScore(b *testing.B) {
	netFile := "../testdata/benchmark/network.nwk"
	geneTrees := "../testdata/benchmark/gene-trees.nex"
	tre, genes, err := prep.ReadInputFiles(netFile, geneTrees, "nexus")
	if err != nil {
		b.Fatalf("failed to read in input files %s", err)
	}
	network, err := prep.ConvertToNetwork(tre)
	if err != nil {
		b.Fatalf("failed to convert tree to network %s", err)
	}
	for i := 0; i < b.N; i++ {
		CalculateReticulationScore(network, genes.Trees)
	}
}
