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
		nQuartets   int
		expectedErr error
	}{
		{
			name:        "basic",
			treeFile:    "../testdata/prep/constraint.nwk",
			quartetFile: "../testdata/prep/quartets.nwk",
			taxaset:     []string{"A", "B", "C", "D", "E", "F", "G", "H", "I", "J"},
			nQuartets:   2,
			expectedErr: nil,
		},
		{
			name:        "non-unique const tree",
			treeFile:    "../testdata/prep/quartets.nwk",
			quartetFile: "../testdata/prep/quartets.nwk",
			taxaset:     []string{},
			nQuartets:   -1,
			expectedErr: ErrInvalidTreeFile,
		},
		{
			name:        "bad const tree",
			treeFile:    "../testdata/prep/badtree.nwk",
			quartetFile: "../testdata/prep/quartets.nwk",
			taxaset:     []string{},
			nQuartets:   -1,
			expectedErr: ErrInvalidNewick,
		},
		{
			name:        "bad const tree (no ;)",
			treeFile:    "../testdata/prep/badtree-nosemi.nwk",
			quartetFile: "../testdata/prep/quartets.nwk",
			taxaset:     []string{},
			nQuartets:   -1,
			expectedErr: ErrInvalidNewick,
		},
		{
			name:        "empty const tree",
			treeFile:    "../testdata/prep/empty.nwk",
			quartetFile: "../testdata/prep/quartets.nwk",
			taxaset:     []string{},
			nQuartets:   -1,
			expectedErr: ErrInvalidTreeFile,
		},
		{
			name:        "bad gene tree trees",
			treeFile:    "../testdata/prep/constraint.nwk",
			quartetFile: "../testdata/prep/badtree.nwk",
			taxaset:     []string{},
			nQuartets:   -1,
			expectedErr: ErrInvalidNewick,
		},
		{
			name:        "empty gene trees",
			treeFile:    "../testdata/prep/constraint.nwk",
			quartetFile: "../testdata/prep/empty.nwk",
			taxaset:     []string{},
			nQuartets:   -1,
			expectedErr: ErrInvalidTreeFile,
		},
	}
	for _, test := range testCases {
		t.Run(test.name, func(t *testing.T) {
			tre, quartets, err := ReadInputFiles(test.treeFile, test.quartetFile)
			if !errors.Is(err, test.expectedErr) {
				t.Errorf("Failed with unexpected error %+v", err)
			} else if errors.Is(err, test.expectedErr) && err != nil {
				t.Logf("%s", err)
			} else if test.expectedErr == nil {
				taxaset := tre.AllTipNames()
				if !reflect.DeepEqual(taxaset, test.taxaset) {
					t.Errorf("Taxaset of tree not equal to expected (%v != %v)", taxaset, test.taxaset)
				}
				if test.nQuartets != len(quartets) {
					t.Errorf("Wrong number of quartets read (%d != %d)", len(quartets), test.nQuartets)
				}
			}
		})
	}
}

func TestReadNetworkFile(t *testing.T) {
	testCases := []struct {
		name             string
		networkFile      string
		expNetwork       string
		expReticulations map[string][2]string
	}{
		{
			name:        "basic test",
			networkFile: "../testdata/prep/net.nwk",
			expNetwork:  "(((9,0),(7,(6,(#H0,8h0u)))),((#H2,(12,((3,(14h2w)#H2),10))h2u),((((5,(#H1,13h1u)),((2h1w)#H1,11))h0w)#H0,(1,4))));",
			expReticulations: map[string][2]string{
				"#H0": {"8h0u", "h0w"},
				"#H1": {"13h1u", "2h1w"},
				"#H2": {"h2u", "14h2w"},
			},
		},
	}
	for _, test := range testCases {
		t.Run(test.name, func(t *testing.T) {
			net, err := ReadNetworkFile(test.networkFile)
			if err != nil {
				t.Fatalf("test returned unexpected err %s", err)
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
					id := nameMap[v[i]]
					if net.Reticulations[k][i] != id {
						t.Errorf("%s mapped to %d when it should have mapped to %d", k, net.Reticulations[k][i], id)
					}
				}
			}
		})
	}
}
