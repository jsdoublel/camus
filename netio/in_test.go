package netio

import (
	"errors"
	"reflect"
	"testing"
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
			treeFile:    "../testdata/netio/constraint.nwk",
			quartetFile: "../testdata/netio/quartets.nwk",
			taxaset:     []string{"A", "B", "C", "D", "E", "F", "G", "H", "I", "J"},
			nQuartets:   2,
			expectedErr: nil,
		},
		{
			name:        "non-unique const tree",
			treeFile:    "../testdata/netio/quartets.nwk",
			quartetFile: "../testdata/netio/quartets.nwk",
			taxaset:     []string{},
			nQuartets:   -1,
			expectedErr: ErrInvalidTreeFile,
		},
		{
			name:        "bad const tree",
			treeFile:    "../testdata/netio/badtree.nwk",
			quartetFile: "../testdata/netio/quartets.nwk",
			taxaset:     []string{},
			nQuartets:   -1,
			expectedErr: ErrInvalidNewick,
		},
		{
			name:        "bad const tree (no ;)",
			treeFile:    "../testdata/netio/badtree-nosemi.nwk",
			quartetFile: "../testdata/netio/quartets.nwk",
			taxaset:     []string{},
			nQuartets:   -1,
			expectedErr: ErrInvalidNewick,
		},
		{
			name:        "empty const tree",
			treeFile:    "../testdata/netio/empty.nwk",
			quartetFile: "../testdata/netio/quartets.nwk",
			taxaset:     []string{},
			nQuartets:   -1,
			expectedErr: ErrInvalidTreeFile,
		},
		{
			name:        "bad gene tree trees",
			treeFile:    "../testdata/netio/constraint.nwk",
			quartetFile: "../testdata/netio/badtree.nwk",
			taxaset:     []string{},
			nQuartets:   -1,
			expectedErr: ErrInvalidNewick,
		},
		{
			name:        "empty gene trees",
			treeFile:    "../testdata/netio/constraint.nwk",
			quartetFile: "../testdata/netio/empty.nwk",
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
			} else if errors.Is(err, test.expectedErr) {
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
