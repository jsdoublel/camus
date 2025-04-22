package prep

import (
	"errors"
	"reflect"
	"strings"
	"testing"

	"github.com/evolbioinfo/gotree/io/newick"
	"github.com/evolbioinfo/gotree/tree"

	gr "github.com/jsdoublel/camus/graphs"
)

func TestIsBinary(t *testing.T) {
	testCases := []struct {
		name     string
		tre      string
		expected bool
	}{
		{
			name:     "basic true",
			tre:      "((a,b),(c,d));",
			expected: true,
		},
		{
			name:     "basic false",
			tre:      "(a,b,c,d);",
			expected: false,
		},
		{
			name:     "unifurcation",
			tre:      "(((a,b)),(c,d));",
			expected: false,
		},
	}
	for _, test := range testCases {
		t.Run(test.name, func(t *testing.T) {
			tre, err := newick.NewParser(strings.NewReader(test.tre)).Parse()
			if err != nil {
				t.Fatal("invalid newick tree; test is written wrong")
			}
			if !TreeIsBinary(tre) == test.expected {
				t.Errorf("got %t expected %t", !TreeIsBinary(tre), test.expected)
			}
		})
	}
}

func TestPreprocess_Errors(t *testing.T) {
	testCases := []struct {
		name        string
		tre         string
		gtrees      []string
		expectedErr error
	}{
		{
			name:        "unrooted",
			tre:         "((a,b),(c,d),e);",
			gtrees:      []string{},
			expectedErr: ErrUnrooted,
		},
		{
			name:        "non-binary",
			tre:         "((a,b),(c,d,e));",
			gtrees:      []string{},
			expectedErr: ErrNonBinary,
		},
		{
			name:        "multree",
			tre:         "((a,a),(a,a));",
			gtrees:      []string{},
			expectedErr: ErrMulTree,
		},
		{
			name:        "multree gtree",
			tre:         "((a,b),(c,d));",
			gtrees:      []string{"((a,a),(a,a));"},
			expectedErr: ErrMulTree,
		},
		{
			name: "missing const labels",
			tre:  "((a,b),(c,d));",
			gtrees: []string{
				"((a,b),(c,d));",
				"(((a,b),(c,d)),e);",
			},
			expectedErr: gr.ErrTipNameMismatch,
		},
		{
			name: "non-binary input trees",
			tre:  "((a,b),(c,d));",
			gtrees: []string{
				"(a,b,c,d);",
			},
			expectedErr: nil,
		},
	}
	for _, test := range testCases {
		t.Run(test.name, func(t *testing.T) {
			tre, err := newick.NewParser(strings.NewReader(test.tre)).Parse()
			if err != nil {
				t.Fatal("invalid newick tree; test is written wrong")
			}
			gtrees := make([]*tree.Tree, len(test.gtrees))
			for i, gt := range test.gtrees {
				tmp, err := newick.NewParser(strings.NewReader(gt)).Parse()
				if err != nil {
					t.Fatal("invalid newick tree; test is written wrong")
				}
				gtrees[i] = tmp
			}
			_, err = Preprocess(tre, gtrees)
			if err != nil && !errors.Is(err, test.expectedErr) {
				t.Errorf("unexpected error %v", err)
			} else if err != nil {
				t.Logf("%s", err)
			}
		})
	}
}

func TestProcessQuartets(t *testing.T) {
	testCases := []struct {
		name     string
		tre      string
		rqList   []string
		expected []string
	}{
		{
			name: "basic",
			tre:  "((((a,b),c),d),f);",
			rqList: []string{
				"(((a,b),c),d);",
				"(((a,b),c),f);",
				"(((a,b),d),f);",
				"(((c,d),f),a);",
				"(((d,b),a),f);",
				"((c,f),(d,b));",
			},
			expected: []string{
				"(((c,d),f),a);",
				"((b,d),(a,f));",
				"((c,f),(d,b));",
			},
		},
	}
	for _, test := range testCases {
		t.Run(test.name, func(t *testing.T) {
			tre, err := newick.NewParser(strings.NewReader(test.tre)).Parse()
			if err != nil {
				t.Fatal("invalid newick tree; test is written wrong")
			}
			tre.UpdateTipIndex()
			rqList := []*tree.Tree{}
			for _, nwk := range test.rqList {
				tr, err := newick.NewParser(strings.NewReader(nwk)).Parse()
				if err != nil {
					t.Errorf("invalid newick tree %s; test is written wrong", nwk)
				}
				rqList = append(rqList, tr)
			}
			result, err := processQuartets(rqList, tre)
			if err != nil {
				t.Errorf("produced error %+v", err)
			}
			expectedList := []*gr.Quartet{}
			for _, nwk := range test.expected {
				tr, err := newick.NewParser(strings.NewReader(nwk)).Parse()
				if err != nil {
					t.Errorf("invalid newick tree %s; test is written wrong", nwk)
				}
				q, err := gr.NewQuartet(tr, tre)
				if err != nil {
					t.Errorf("invalid newick tree %s; test is written wrong", nwk)
				}
				expectedList = append(expectedList, q)
			}
			expected := make(map[gr.Quartet]uint)
			for _, q := range expectedList {
				expected[*q] += 1
			}
			if !reflect.DeepEqual(*result, expected) {
				t.Errorf("actual %s != expected %s", gr.QSetToString(*result, tre), gr.QSetToString(expected, tre))
			}
		})
	}
}
