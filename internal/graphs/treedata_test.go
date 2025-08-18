package graphs

import (
	"fmt"
	"strings"
	"testing"

	"github.com/bits-and-blooms/bitset"
	"github.com/evolbioinfo/gotree/io/newick"
	"github.com/evolbioinfo/gotree/tree"
)

func TestMakeTreeData(t *testing.T) {
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
			err = tre.UpdateTipIndex()
			if err != nil {
				t.Error(err)
			}
			qc, err := makeQCounts(q, tre)
			if err != nil {
				t.Error("invalid quartet; test is written wrong")
			}

			treeData := MakeTreeData(tre, qc)
			lca := treeData.lca
			leafset := treeData.leafsets
			quartetSets := treeData.quartetSet
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

func quartetSetEqualityTester(quartetSets [][]Quartet, testQS map[string][]string, tre *tree.Tree) (bool, error) {
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
			q1, err := NewQuartet(qTree, tre)
			if err != nil {
				return false, err
			}
			found := false
			for _, q2 := range quartetSets[node.Id()] {
				if q1.Compare(q2) == Qeq {
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

func makeQCounts(qList []*tree.Tree, constTree *tree.Tree) (map[Quartet]uint, error) {
	result := make(map[Quartet]uint)
	for _, qt := range qList {
		q, err := NewQuartet(qt, constTree)
		if err != nil {
			return nil, err
		}
		result[q] += 1
	}
	return result, nil
}
