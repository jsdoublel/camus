package prep

import (
	"errors"
	"fmt"
	"reflect"
	"strings"
	"testing"

	"github.com/bits-and-blooms/bitset"
	"github.com/evolbioinfo/gotree/io/newick"
	"github.com/evolbioinfo/gotree/tree"

	"github.com/jsdoublel/camus/qrt"
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
			expectedErr: ErrInvalidTree,
		},
		{
			name:        "non-binary",
			tre:         "(a,b,c,d);",
			gtrees:      []string{},
			expectedErr: ErrInvalidTree,
		},
		{
			name:        "multree",
			tre:         "((a,a),(a,a));",
			gtrees:      []string{},
			expectedErr: ErrInvalidTree,
		},
		{
			name:        "multree gtree",
			tre:         "((a,b),(c,d));",
			gtrees:      []string{"((a,a),(a,a));"},
			expectedErr: ErrInvalidTree,
		},
		{
			name: "missing const labels",
			tre:  "((a,b),(c,d));",
			gtrees: []string{
				"((a,b),(c,d));",
				"(((a,b),(c,d)),e);",
			},
			expectedErr: qrt.ErrTipNameMismatch,
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
			_, _, err = Preprocess(tre, gtrees)
			if !errors.Is(err, test.expectedErr) {
				t.Errorf("unexpected error %v", err)
			} else {
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
			result, _, err := processQuartets(rqList, tre)
			if err != nil {
				t.Errorf("produced error %+v", err)
			}
			expectedList := []*qrt.Quartet{}
			for _, nwk := range test.expected {
				tr, err := newick.NewParser(strings.NewReader(nwk)).Parse()
				if err != nil {
					t.Errorf("invalid newick tree %s; test is written wrong", nwk)
				}
				q, err := qrt.NewQuartet(tr, tre)
				if err != nil {
					t.Errorf("invalid newick tree %s; test is written wrong", nwk)
				}
				expectedList = append(expectedList, q)
			}
			expected := make(map[qrt.Quartet]uint)
			for _, q := range expectedList {
				expected[*q] += 1
			}
			if !reflect.DeepEqual(*result, expected) {
				t.Errorf("actual %s != expected %s", qrt.QSetToString(*result, tre), qrt.QSetToString(expected, tre))
			}
		})
	}
}

func TestProcessTreeData(t *testing.T) {
	testCases := []struct {
		name        string
		tre         string
		quartets    []string
		lca         map[string][][]string
		leafset     map[string][]string
		quartetSets map[string][]string
	}{
		{
			name:     "basic",
			tre:      "((((A,B)a,C)b,D)c,F)r;",
			quartets: []string{"((A,C),(B,D));"},
			lca: map[string][][]string{
				"a": {
					{"A", "B"},
				},
				"b": {
					{"A", "C"},
					{"B", "C"},
				},
				"c": {
					{"A", "D"},
					{"B", "D"},
					{"C", "D"},
				},
				"r": {
					{"A", "F"},
					{"B", "F"},
					{"C", "F"},
					{"D", "F"},
				},
			},
			leafset: map[string][]string{
				"a": {"A", "B"},
				"b": {"A", "B", "C"},
				"c": {"A", "B", "C", "D"},
				"r": {"A", "B", "C", "D", "F"},
			},
			quartetSets: map[string][]string{
				"a": {},
				"b": {},
				"c": {"((A,C),(B,D));"},
				"r": {"((A,C),(B,D));"},
			},
		},
	}
	for _, test := range testCases {
		t.Run(test.name, func(t *testing.T) {
			tre, err := newick.NewParser(strings.NewReader(test.tre)).Parse()
			if err != nil {
				t.Error("invalid newick tree; test is written wrong")
			}
			q := make([]*tree.Tree, 0)
			for _, s := range test.quartets {
				tmp, err := newick.NewParser(strings.NewReader(s)).Parse()
				if err != nil {
					t.Error("invalid newick tree; test is written wrong")
				}
				q = append(q, tmp)
			}
			tre.UpdateTipIndex()
			qs, _, err := processQuartets(q, tre)
			treeData := MakeTreeData(tre, qs)
			lca := treeData.lca
			leafset := treeData.leafsets
			quartetSets := treeData.QuartetSet
			nLeaves := len(lca)
			for i := range nLeaves {
				for j := range nLeaves {
					if lca[i][j] != lca[j][i] {
						t.Error("lca structure problem")
					}
				}
			}
			if b, err := lcaEqualityTester(lca, test.lca, tre); err != nil {
				t.Error(err.Error())
			} else if !b {
				t.Error("lca != expected")
			}
			if b, err := leafsetEqualityTester(leafset, test.leafset, tre); err != nil {
				t.Error(err.Error())
			} else if !b {
				t.Error("leafset != expected")
			}
			if b, err := quartetSetEqualityTester(quartetSets, test.quartetSets, tre); err != nil {
				t.Error(err.Error())
			} else if !b {
				t.Error("quartetSets != expected")
			}

		})
	}
}

func lcaEqualityTester(lca [][]int, testLCA map[string][][]string, tre *tree.Tree) (bool, error) {
	for k, v := range testLCA {
		for _, pair := range v {
			if len(pair) != 2 {
				return false, fmt.Errorf("lca is not a pair; test is written wrong")
			}
			node1, err := getNode(pair[0], tre)
			if err != nil {
				return false, err
			}
			node2, err := getNode(pair[1], tre)
			if err != nil {
				return false, err
			}
			lcaNode, err := getNode(k, tre)
			if err != nil {
				return false, err
			}
			if lca[node1.Id()][node2.Id()] != lcaNode.Id() {
				return false, nil
			}
		}
	}
	return true, nil
}

func leafsetEqualityTester(leafset []*bitset.BitSet, testLeafset map[string][]string, tre *tree.Tree) (bool, error) {
	for k, v := range testLeafset {
		node, err := getNode(k, tre)
		if err != nil {
			return false, err
		}
		leafsetList := leafset[node.Id()]
		for _, leaf := range v {
			id, err := tre.TipIndex(leaf)
			if err != nil {
				return false, err
			}
			if !leafsetList.Test(uint(id)) {
				return false, nil
			}
		}
	}
	return true, nil
}

func quartetSetEqualityTester(quartetSets [][]*qrt.Quartet, testQS map[string][]string, tre *tree.Tree) (bool, error) {
	for k, v := range testQS {
		node, err := getNode(k, tre)
		if err != nil {
			return false, err
		}
		for _, quartetString := range v {
			qTree, err := newick.NewParser(strings.NewReader(quartetString)).Parse()
			if err != nil {
				return false, fmt.Errorf("cannot parse %s as newick tree. %w", quartetString, err)
			}
			q1, err := qrt.NewQuartet(qTree, tre)
			if err != nil {
				return false, err
			}
			found := false
			for _, q2 := range quartetSets[node.Id()] {
				if q1.Compare(q2) == qrt.Qeq {
					found = true
				}
			}
			if !found {
				return false, nil
			}
		}
	}
	return true, nil
}

func getNode(label string, tre *tree.Tree) (*tree.Node, error) {
	nodeList, err := tre.SelectNodes(label)
	if err != nil {
		panic(err)
	}
	if len(nodeList) != 1 {
		return nil, fmt.Errorf("more or less than one internal node with the required label; test is written wrong")
	}
	return nodeList[0], nil
}
