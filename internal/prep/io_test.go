package prep

import (
	"errors"
	"path/filepath"
	"reflect"
	"testing"

	"github.com/evolbioinfo/gotree/tree"
)

const testData = "../../testdata"

func TestReadInputFiles(t *testing.T) {
	testCases := []struct {
		name        string
		treeFile    string
		quartetFile string
		taxaset     []string
		numGenes    int
		format      string
		expectedErr error
	}{
		{
			name:        "basic",
			treeFile:    filepath.Join(testData, "prep/constraint.nwk"),
			quartetFile: filepath.Join(testData, "prep/quartets.nwk"),
			taxaset:     []string{"A", "B", "C", "D", "E", "F", "G", "H", "I", "J"},
			numGenes:    2,
			format:      "newick",
			expectedErr: nil,
		},
		{
			name:        "non-unique const tree",
			treeFile:    filepath.Join(testData, "prep/quartets.nwk"),
			quartetFile: filepath.Join(testData, "prep/quartets.nwk"),
			taxaset:     []string{},
			numGenes:    -1,
			format:      "newick",
			expectedErr: ErrInvalidFile,
		},
		{
			name:        "bad const tree",
			treeFile:    filepath.Join(testData, "prep/badtree.nwk"),
			quartetFile: filepath.Join(testData, "prep/quartets.nwk"),
			taxaset:     []string{},
			numGenes:    -1,
			format:      "newick",
			expectedErr: ErrInvalidFormat,
		},
		{
			name:        "bad const tree (no ;)",
			treeFile:    filepath.Join(testData, "prep/badtree-nosemi.nwk"),
			quartetFile: filepath.Join(testData, "prep/quartets.nwk"),
			taxaset:     []string{},
			numGenes:    -1,
			format:      "newick",
			expectedErr: ErrInvalidFormat,
		},
		{
			name:        "empty const tree",
			treeFile:    filepath.Join(testData, "prep/empty.nwk"),
			quartetFile: filepath.Join(testData, "prep/quartets.nwk"),
			taxaset:     []string{},
			numGenes:    -1,
			format:      "newick",
			expectedErr: ErrInvalidFile,
		},
		{
			name:        "bad gene tree trees",
			treeFile:    filepath.Join(testData, "prep/constraint.nwk"),
			quartetFile: filepath.Join(testData, "prep/badtree.nwk"),
			taxaset:     []string{},
			numGenes:    -1,
			format:      "newick",
			expectedErr: ErrInvalidFormat,
		},
		{
			name:        "empty gene trees",
			treeFile:    filepath.Join(testData, "prep/constraint.nwk"),
			quartetFile: filepath.Join(testData, "prep/empty.nwk"),
			taxaset:     []string{},
			numGenes:    -1,
			format:      "newick",
			expectedErr: ErrInvalidFile,
		},
		{
			name:        "basic nexus",
			treeFile:    filepath.Join(testData, "prep/constraint.nwk"),
			quartetFile: filepath.Join(testData, "prep/quartets.nex"),
			taxaset:     []string{"A", "B", "C", "D", "E", "F", "G", "H", "I", "J"},
			numGenes:    2,
			format:      "nexus",
			expectedErr: nil,
		},
	}
	for _, test := range testCases {
		t.Run(test.name, func(t *testing.T) {
			tre, quartets, err := ReadInputFiles(test.treeFile, test.quartetFile, parseFormat[test.format])
			switch {
			case !errors.Is(err, test.expectedErr):
				t.Errorf("Failed with unexpected error %+v", err)
			case errors.Is(err, test.expectedErr) && err != nil:
				t.Logf("%s", err)
			case test.expectedErr == nil:
				taxaset := tre.AllTipNames()
				if !reflect.DeepEqual(taxaset, test.taxaset) {
					t.Errorf("Taxaset of tree not equal to expected (%v != %v)", taxaset, test.taxaset)
				}
				if test.numGenes != len(quartets.Trees) {
					t.Errorf("Wrong number of quartets read (%d != %d)", len(quartets.Trees), test.numGenes)
				}
			}
		})
	}
}

func TestConvertToNetwork(t *testing.T) {
	testCases := []struct {
		name             string
		networkFile      string
		expNetwork       string
		expReticulations map[string][2]string
		expectedErr      error
	}{
		{
			name:        "basic test",
			networkFile: filepath.Join(testData, "prep/net.nwk"),
			expNetwork:  "(((9,0),(7,(6,(#H0,8h0u)))),((#H2,(12,((3,(14h2w)#H2),10))h2u),((((5,(#H1,13h1u)),((2h1w)#H1,11))h0w)#H0,(1,4))));",
			expReticulations: map[string][2]string{
				"#H0": {"8h0u", "h0w"},
				"#H1": {"13h1u", "2h1w"},
				"#H2": {"h2u", "14h2w"},
			},
			expectedErr: nil,
		},
		{
			name:             "unresolved",
			networkFile:      filepath.Join(testData, "prep/unresolved.nwk"),
			expNetwork:       "",
			expReticulations: nil,
			expectedErr:      ErrNonBinary,
		},
		{
			name:             "non-unique network",
			networkFile:      filepath.Join(testData, "prep/multi-net.nwk"),
			expReticulations: nil,
			expectedErr:      ErrInvalidFile,
		},
		{
			name:             "bad network",
			networkFile:      filepath.Join(testData, "prep/badtree.nwk"),
			expReticulations: nil,
			expectedErr:      ErrInvalidFormat,
		},
		{
			name:             "bad network (no ;)",
			networkFile:      filepath.Join(testData, "prep/badtree-nosemi.nwk"),
			expReticulations: nil,
			expectedErr:      ErrInvalidFormat,
		},
		{
			name:             "empty",
			networkFile:      filepath.Join(testData, "prep/empty.nwk"),
			expReticulations: nil,
			expectedErr:      ErrInvalidFile,
		},
		{
			name:             "no reticulations",
			networkFile:      filepath.Join(testData, "prep/constraint.nwk"),
			expReticulations: nil,
			expectedErr:      ErrNoReticulations,
		},
		{
			name:             "unrooted",
			networkFile:      filepath.Join(testData, "prep/unrooted-net.nwk"),
			expReticulations: nil,
			expectedErr:      ErrUnrooted,
		},
	}
	for _, test := range testCases {
		t.Run(test.name, func(t *testing.T) {
			tre, err := readTreeFile(test.networkFile)
			if err != nil && !errors.Is(err, test.expectedErr) {
				t.Fatalf("test returned unexpected err %s", err)
			} else if err != nil && errors.Is(err, test.expectedErr) {
				t.Logf("%s", err)
				return
			}
			net, err := ConvertToNetwork(tre)
			if err != nil && !errors.Is(err, test.expectedErr) {
				t.Fatalf("test returned unexpected err %s", err)
			} else if err != nil && errors.Is(err, test.expectedErr) {
				t.Logf("%s", err)
				return
			}
			if net.Newick() != test.expNetwork {
				t.Errorf("result != expected, %s != %s", net.Newick(), test.expNetwork)
			}
			nameMap := make(map[string]int)
			net.NetTree.PostOrder(func(cur, prev *tree.Node, e *tree.Edge) (keep bool) {
				nameMap[cur.Name()] = cur.Id()
				return true
			})
			for k, v := range test.expReticulations {
				for i := range 2 {
					if id := nameMap[v[i]]; net.Reticulations[k].IDs[i] != id {
						t.Errorf("%s mapped to %d when it should have mapped to %d", k, net.Reticulations[k].IDs[i], id)
					}
				}
			}
		})
	}
}
