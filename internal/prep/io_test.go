package prep

import (
	"errors"
	"reflect"
	"testing"

	"github.com/evolbioinfo/gotree/tree"
)

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
			treeFile:    "testdata/constraint.nwk",
			quartetFile: "testdata/quartets.nwk",
			taxaset:     []string{"A", "B", "C", "D", "E", "F", "G", "H", "I", "J"},
			numGenes:    2,
			format:      "newick",
			expectedErr: nil,
		},
		{
			name:        "non-unique const tree",
			treeFile:    "testdata/quartets.nwk",
			quartetFile: "testdata/quartets.nwk",
			taxaset:     []string{},
			numGenes:    -1,
			format:      "newick",
			expectedErr: ErrInvalidFile,
		},
		{
			name:        "bad const tree",
			treeFile:    "testdata/badtree.nwk",
			quartetFile: "testdata/quartets.nwk",
			taxaset:     []string{},
			numGenes:    -1,
			format:      "newick",
			expectedErr: ErrInvalidFormat,
		},
		{
			name:        "bad const tree (no ;)",
			treeFile:    "testdata/badtree-nosemi.nwk",
			quartetFile: "testdata/quartets.nwk",
			taxaset:     []string{},
			numGenes:    -1,
			format:      "newick",
			expectedErr: ErrInvalidFormat,
		},
		{
			name:        "empty const tree",
			treeFile:    "testdata/empty.nwk",
			quartetFile: "testdata/quartets.nwk",
			taxaset:     []string{},
			numGenes:    -1,
			format:      "newick",
			expectedErr: ErrInvalidFile,
		},
		{
			name:        "bad gene tree trees",
			treeFile:    "testdata/constraint.nwk",
			quartetFile: "testdata/badtree.nwk",
			taxaset:     []string{},
			numGenes:    -1,
			format:      "newick",
			expectedErr: ErrInvalidFormat,
		},
		{
			name:        "empty gene trees",
			treeFile:    "testdata/constraint.nwk",
			quartetFile: "testdata/empty.nwk",
			taxaset:     []string{},
			numGenes:    -1,
			format:      "newick",
			expectedErr: ErrInvalidFile,
		},
		{
			name:        "basic nexus",
			treeFile:    "testdata/constraint.nwk",
			quartetFile: "testdata/quartets.nex",
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
			networkFile: "testdata/net.nwk",
			expNetwork:  "(((9,0),(7,(6,(#H1,8h0u)))),((#H3,(12,((3,(14h2w)#H3),10))h2u),((((5,(#H2,13h1u)),((2h1w)#H2,11))h0w)#H1,(1,4))));",
			expReticulations: map[string][2]string{
				"#H1": {"8h0u", "h0w"},
				"#H2": {"13h1u", "2h1w"},
				"#H3": {"h2u", "14h2w"},
			},
			expectedErr: nil,
		},
		{
			name:             "unresolved",
			networkFile:      "testdata/unresolved.nwk",
			expNetwork:       "",
			expReticulations: nil,
			expectedErr:      ErrNonBinary,
		},
		{
			name:             "non-unique network",
			networkFile:      "testdata/multi-net.nwk",
			expReticulations: nil,
			expectedErr:      ErrInvalidFile,
		},
		{
			name:             "bad network",
			networkFile:      "testdata/badtree.nwk",
			expReticulations: nil,
			expectedErr:      ErrInvalidFormat,
		},
		{
			name:             "bad network (no ;)",
			networkFile:      "testdata/badtree-nosemi.nwk",
			expReticulations: nil,
			expectedErr:      ErrInvalidFormat,
		},
		{
			name:             "empty",
			networkFile:      "testdata/empty.nwk",
			expReticulations: nil,
			expectedErr:      ErrInvalidFile,
		},
		{
			name:             "no reticulations",
			networkFile:      "testdata/constraint.nwk",
			expReticulations: nil,
			expectedErr:      ErrNoReticulations,
		},
		{
			name:             "unrooted",
			networkFile:      "testdata/unrooted-net.nwk",
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
