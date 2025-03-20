package netio

import (
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
	}{
		{
			name:        "basic",
			treeFile:    "../testdata/netio/constraint.nwk",
			quartetFile: "../testdata/netio/quartets.nwk",
			taxaset:     []string{"A", "B", "C", "D", "E", "F", "G", "H", "I", "J"},
			nQuartets:   2,
		},
	}
	for _, test := range testCases {
		t.Run(test.name, func(t *testing.T) {
			tre, quartets, err := ReadInputFiles(test.treeFile, test.quartetFile)
			if err != nil {
				t.Errorf("Failed with error %+v", err)
			}
			taxaset := tre.AllTipNames()
			if !reflect.DeepEqual(taxaset, test.taxaset) {
				t.Errorf("Taxaset of tree not equal to expected (%v != %v)", taxaset, test.taxaset)
			}
			if test.nQuartets != len(quartets) {
				t.Errorf("Wrong number of quartets read (%d != %d)", len(quartets), test.nQuartets)
			}
		})
	}
}
