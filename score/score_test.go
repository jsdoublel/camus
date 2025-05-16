package score

import (
	"bytes"
	"errors"
	"io"
	"math"
	"os"
	"strings"
	"testing"

	"github.com/evolbioinfo/gotree/io/newick"
	"github.com/evolbioinfo/gotree/tree"

	pr "github.com/jsdoublel/camus/prep"
)

func TestCalculateRecticulationScore(t *testing.T) {
	testCases := []struct {
		name        string
		network     string
		gtrees      []string
		expected    []*map[string]float64
		expectedErr error
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
			expectedErr: nil,
		},
		{
			name:        "not level-1",
			network:     "(A,(B,(#H1,(C,(#H0,(D,(E,(F,((G,(H,((I,J))#H1)))#H0))))))));",
			gtrees:      nil,
			expected:    nil,
			expectedErr: ErrNotLevel1,
		},
	}
	for _, test := range testCases {
		t.Run(test.name, func(t *testing.T) {
			tre, err := newick.NewParser(strings.NewReader(test.network)).Parse()
			if err != nil {
				t.Fatalf("invalid newick in file %s", err)
			}
			ntw, err := pr.ConvertToNetwork(tre)
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
			switch {
			case err != nil && !errors.Is(err, test.expectedErr):
				t.Errorf("test case failed with unexpected error %s", err)
			case err != nil:
				t.Logf("%s", err)
			case !compareScoreMaps(result, test.expected):
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

func TestCalculateRecticulationScore_Large(t *testing.T) {
	testCases := []struct {
		name     string
		network  string
		gtrees   string
		expected string
	}{
		{
			name:     "pauls data",
			network:  "../testdata/large/network.nwk",
			gtrees:   "../testdata/large/gene-trees.nwk",
			expected: "../testdata/large/scores.csv",
		},
	}
	for _, test := range testCases {
		t.Run(test.name, func(t *testing.T) {
			tre, genes, err := pr.ReadInputFiles(test.network, test.gtrees, pr.Newick)
			if err != nil {
				t.Fatalf("failed to read in input files %s", err)
			}
			network, err := pr.ConvertToNetwork(tre)
			if err != nil {
				t.Fatalf("failed to convert tree to network %s", err)
			}
			scores, err := CalculateReticulationScore(network, genes.Trees)
			if err != nil {
				t.Fatalf("failed with unexpected err %s", err)
			}
			r, w, err := os.Pipe()
			if err != nil {
				t.Fatalf("failed to open pipe %s", err)
			}
			oldStdout := os.Stdout
			os.Stdout = w
			pr.WriteBranchScoresToCSV(scores, genes.Names)
			err = w.Close()
			if err != nil {
				t.Fatalf("could not close pipe: %s", err)
			}
			os.Stdout = oldStdout
			var buf bytes.Buffer
			_, err = io.Copy(&buf, r)
			if err != nil {
				t.Fatalf("could not copy stdout: %s", err)
			}
			result := strings.TrimSpace(buf.String())
			expBytes, err := os.ReadFile(test.expected)
			if err != nil {
				t.Fatalf("failed to read in expected file %s", err)
			}
			if result != strings.TrimSpace(string(expBytes)) {
				t.Error("result != expected")
			}
		})
	}
}

func BenchmarkCalculateRecticulationScore(b *testing.B) {
	netFile := "../testdata/large/network.nwk"
	geneTrees := "../testdata/large/gene-trees.nwk"
	tre, genes, err := pr.ReadInputFiles(netFile, geneTrees, pr.Newick)
	if err != nil {
		b.Fatalf("failed to read in input files %s", err)
	}
	network, err := pr.ConvertToNetwork(tre)
	if err != nil {
		b.Fatalf("failed to convert tree to network %s", err)
	}
	for b.Loop() {
		_, err := CalculateReticulationScore(network, genes.Trees)
		if err != nil {
			b.Fatalf("Failed to calculate reticulation scores: %s", err)
		}
	}
}
