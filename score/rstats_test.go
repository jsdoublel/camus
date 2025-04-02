package score

import (
	"reflect"
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
			network: "../testdata/prep/net.nwk",
			gtrees: []string{
				"((0,9),(5,7));",
				"((5,7),(9,6));",
				"((9,7),(5,6));",
				"((5,9),(7,6));",
			},
			expected: []*map[string]float64{
				{"#H0": float64(0), "#H1": float64(0), "#H2": float64(0)},
				{"#H0": float64(0), "#H1": float64(0), "#H2": float64(0)},
				{"#H0": float64(1), "#H1": float64(0), "#H2": float64(0)},
				{"#H0": float64(0), "#H1": float64(0), "#H2": float64(0)},
			},
		},
	}
	for _, test := range testCases {
		t.Run(test.name, func(t *testing.T) {
			ntw, err := prep.ReadNetworkFile(test.network)
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
			if !reflect.DeepEqual(result, test.expected) {
				t.Errorf("%v != %v; result != expected", result, test.expected)
			}

		})
	}
}
