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

	pr "github.com/jsdoublel/camus/internal/prep"
)

func TestRecticulationScore(t *testing.T) {
	testCases := []struct {
		name        string
		network     string
		gtrees      []string
		expected    []*map[string]float64
		expectedErr error
	}{
		{
			name:    "basic",
			network: "(((9,0),(7,(6,(#H1,8h0u)))),((#H3,(12,((3,(14h2w)#H3),10))h2u),((((5,(#H2,13h1u)),((2h1w)#H2,11))h0w)#H1,(1,4))));",
			gtrees: []string{
				"((0,9),(5,7));",
				"((5,7),(9,6));",
				"((9,7),(5,6));",
				"((5,9),(7,6));",
			},
			expected: []*map[string]float64{
				{"#H1": math.NaN(), "#H2": math.NaN(), "#H3": math.NaN()},
				{"#H1": float64(0), "#H2": math.NaN(), "#H3": math.NaN()},
				{"#H1": float64(1), "#H2": math.NaN(), "#H3": math.NaN()},
				{"#H1": float64(0), "#H2": math.NaN(), "#H3": math.NaN()},
			},
			expectedErr: nil,
		},
		{
			name:        "not level-1",
			network:     "(A,(B,(#H2,(C,(#H1,(D,(E,(F,((G,(H,((I,J))#H2)))#H1))))))));",
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
			result, err := ReticulationScore(ntw, gtrees)
			switch {
			case err != nil && !errors.Is(err, test.expectedErr):
				t.Errorf("test case failed with unexpected error %s", err)
			case err != nil:
				t.Logf("%s", err)
			default:
				compareScoreMaps(t, result, test.expected)
			}
		})
	}
}

// compares the two maps (specifically allows NaN == NaN to be true)
func compareScoreMaps(t *testing.T, got, want []*map[string]float64) {
	t.Helper()
	if len(got) != len(want) {
		t.Fatalf("score map length %d != %d", len(got), len(want))
	}
	for i := range got {
		if len(*got[i]) != len(*want[i]) {
			t.Fatalf("entry %d map size %d != %d", i, len(*got[i]), len(*want[i]))
		}
		for key := range *got[i] {
			wantVal, ok := (*want[i])[key]
			if !ok {
				t.Fatalf("missing key %q in expected map index %d", key, i)
			}
			gotVal := (*got[i])[key]
			if gotVal != wantVal && (!math.IsNaN(gotVal) || !math.IsNaN(wantVal)) {
				t.Fatalf("map[%d][%q] = %v, want %v", i, key, gotVal, wantVal)
			}
		}
	}
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
			network:  "testdata/network.nwk",
			gtrees:   "testdata/gene-trees.nwk",
			expected: "testdata/scores.csv",
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
			scores, err := ReticulationScore(network, genes.Trees)
			if err != nil {
				t.Fatalf("failed with unexpected err %s", err)
			}
			r, w, err := os.Pipe()
			if err != nil {
				t.Fatalf("failed to open pipe %s", err)
			}
			oldStdout := os.Stdout
			os.Stdout = w
			pr.WriteRetScoresToCSV(scores, genes.Names)
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
	netFile := "testdata/network.nwk"
	geneTrees := "testdata/gene-trees.nwk"
	tre, genes, err := pr.ReadInputFiles(netFile, geneTrees, pr.Newick)
	if err != nil {
		b.Fatalf("failed to read in input files %s", err)
	}
	network, err := pr.ConvertToNetwork(tre)
	if err != nil {
		b.Fatalf("failed to convert tree to network %s", err)
	}
	for b.Loop() {
		_, err := ReticulationScore(network, genes.Trees)
		if err != nil {
			b.Fatalf("Failed to calculate reticulation scores: %s", err)
		}
	}
}
