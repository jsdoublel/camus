package prep

import (
	"fmt"
	"reflect"
	"slices"
	"strings"
	"testing"

	"github.com/evolbioinfo/gotree/io/newick"
	"github.com/evolbioinfo/gotree/tree"
)

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
				t.Error("invalid newick tree; test is written wrong")
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
			expected := []*Quartet{}
			for _, nwk := range test.expected {
				tr, err := newick.NewParser(strings.NewReader(nwk)).Parse()
				if err != nil {
					t.Errorf("invalid newick tree %s; test is written wrong", nwk)
				}
				q, err := NewQuartet(tr, tre)
				if err != nil {
					t.Errorf("invalid newick tree %s; test is written wrong", nwk)
				}
				expected = append(expected, q)
			}
			slices.SortFunc(result, sortQuartet)
			slices.SortFunc(expected, sortQuartet)
			if !reflect.DeepEqual(result, expected) {
				t.Errorf("actual %s != expected %s", listToString(result, tre), listToString(expected, tre))
			}
		})
	}
}

func sortQuartet(q1, q2 *Quartet) int {
	sum1 := 0
	sum2 := 0
	for i := 0; i < 4; i++ {
		sum1 += int(q1.taxa[i])
		sum2 += int(q2.taxa[i])
	}
	sum1 += int(q1.topology)
	sum2 += int(q2.topology)
	return sum1 - sum2
}

func TestProcessTreeData(t *testing.T) {
	testCases := []struct {
		name    string
		tre     string
		lca     map[string][][]string
		leafset map[string][]string
	}{
		{
			name: "basic",
			tre:  "((((A,B)a,C)b,D)c,F)r;",
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
		},
	}
	for _, test := range testCases {
		t.Run(test.name, func(t *testing.T) {
			tre, err := newick.NewParser(strings.NewReader(test.tre)).Parse()
			if err != nil {
				t.Error("invalid newick tree; test is written wrong")
			}
			tre.UpdateTipIndex()
			treeData := PreprocessTreeData(tre)
			lca := treeData.LCA
			leafset := treeData.Leafsets
			nLeaves := len(lca)
			for i := 0; i < nLeaves; i++ {
				for j := 0; j < nLeaves; j++ {
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
		})
	}
}

func lcaEqualityTester(lca [][]uint, testLCA map[string][][]string, tre *tree.Tree) (bool, error) {
	for k, v := range testLCA {
		for _, pair := range v {
			if len(pair) != 2 {
				return false, fmt.Errorf("lca is not a pair; test is written wrong")
			}
			id1, err := tre.TipIndex(pair[0])
			tipIndexPanic("lca test", err)
			id2, err := tre.TipIndex(pair[1])
			tipIndexPanic("lca test", err)
			nodeList, err := tre.SelectNodes(k)
			if err != nil {
				panic(err)
			}
			if len(nodeList) != 1 {
				return false, fmt.Errorf("more or less than one internal node with the required label; test is written wrong")
			}
			if lca[id1][id2] != uint(nodeList[0].Id()) {
				return false, nil
			}
		}
	}
	return true, nil
}

func leafsetEqualityTester(leafset [][]bool, testLeafset map[string][]string, tre *tree.Tree) (bool, error) {
	for k, v := range testLeafset {
		nodeList, err := tre.SelectNodes(k)
		if err != nil {
			panic(err)
		}
		if len(nodeList) != 1 {
			return false, fmt.Errorf("more or less than one internal node with the required label; test is written wrong")
		}
		leafsetList := leafset[nodeList[0].Id()]
		for _, leaf := range v {
			id, err := tre.TipIndex(leaf)
			tipIndexPanic("leafset test", err)
			if !leafsetList[id] {
				return false, nil
			}
		}
	}
	return true, nil
}
