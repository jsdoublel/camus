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

type TestQuartet struct {
	set1 []string
	set2 []string
}

func TestNewQuartet(t *testing.T) {
	testCases := []struct {
		name    string
		qTree   string
		tre     string
		quartet *TestQuartet
	}{
		{
			name:  "basic",
			qTree: "((a,b),(c,d));",
			tre:   "(((a,c),(b,d)),f);",
			quartet: &TestQuartet{
				set1: []string{"a", "b"},
				set2: []string{"c", "d"},
			},
		},
		{
			name:  "unbalanced",
			qTree: "(((a, b), c), d);",
			tre:   "(((a,c),(b,d)),f);",
			quartet: &TestQuartet{
				set1: []string{"a", "b"},
				set2: []string{"c", "d"},
			},
		},
	}

	for _, test := range testCases {
		t.Run(test.name, func(t *testing.T) {
			qTree, err := newick.NewParser(strings.NewReader(test.qTree)).Parse()
			if err != nil {
				t.Error("invalid newick tree; test is written wrong")
			}
			tre, err := newick.NewParser(strings.NewReader(test.tre)).Parse()
			tre.UpdateTipIndex()
			if err != nil {
				t.Error("invalid newick tree; test is written wrong")
			}
			q, err := NewQuartet(qTree, tre)
			if err != nil {
				t.Errorf("produced err %+v", err)
			}
			t.Logf("\nquartet %s\ntest quartet %s", q.String(tre), test.quartet.String())
			if !testQuartetEqual(q, test.quartet, tre) {
				t.Errorf("quartet has wrong topology (%s != %s)", q.String(tre), test.quartet.String())
			}
		})
	}
}

func TestCompare(t *testing.T) {
	testCases := []struct {
		name   string
		q1     string
		q2     string
		tre    string
		result int
	}{
		{
			name:   "equals",
			q1:     "((a,b),(c,d));",
			q2:     "(((a, b), c), d);",
			tre:    "(((a,c),(b,d)),f);",
			result: Q_EQ,
		},
		{
			name:   "not equal",
			q1:     "((a,b),(c,d));",
			q2:     "(((a, c), b), d);",
			tre:    "(((a,c),(b,d)),f);",
			result: Q_NEQ,
		},
		{
			name:   "diff",
			q1:     "((a,b),(c,d));",
			q2:     "(((a, f), b), d);",
			tre:    "(((a,c),(b,d)),f);",
			result: Q_DIFF,
		},
	}
	for _, test := range testCases {
		t.Run(test.name, func(t *testing.T) {
			t1, err := newick.NewParser(strings.NewReader(test.q1)).Parse()
			if err != nil {
				t.Error("invalid newick tree; test is written wrong")
			}
			t2, err := newick.NewParser(strings.NewReader(test.q2)).Parse()
			if err != nil {
				t.Error("invalid newick tree; test is written wrong")
			}
			tre, err := newick.NewParser(strings.NewReader(test.tre)).Parse()
			if err != nil {
				t.Error("invalid newick tree; test is written wrong")
			}
			tre.UpdateTipIndex()
			q1, err := NewQuartet(t1, tre)
			if err != nil {
				t.Errorf("produced err %+v", err)
			}
			q2, err := NewQuartet(t2, tre)
			if err != nil {
				t.Errorf("produced err %+v", err)
			}
			t.Logf("\nquartet %s\ntest quartet %s", q1.String(tre), q2.String(tre))
			result := q1.Compare(q2)
			if result != test.result {
				t.Errorf("got %d, expected %d", result, test.result)
			}
		})
	}
}

func TestQuartetsFromTree(t *testing.T) {
	testCases := []struct {
		name string
		tre  string
		qSet []string
	}{
		{
			name: "basic",
			tre:  "((((a,b),c),d),f);",
			qSet: []string{
				"(((a,b),c),d);",
				"(((a,b),c),f);",
				"(((a,b),d),f);",
				"(((a,c),d),f);",
				"(((b,c),d),f);",
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
			qSet := QuartetsFromTree(tre)
			expectedQSet, err := stringListToQMap(test.qSet, tre)
			if err != nil {
				t.Error(err)
			}
			if !reflect.DeepEqual(qSet, expectedQSet) {
				t.Errorf("actual %s != expected %s", setToString(qSet, tre), setToString(expectedQSet, tre))
			}
		})
	}
}

func (tq *TestQuartet) topology(tre *tree.Tree) uint8 {
	ids := make([]int, 4)
	partition := make(map[int]bool)
	for i := 0; i < 4; i++ {
		if i < 2 {
			ti, err := tre.TipIndex(tq.set1[i])
			tipIndexPanic("testQuartet topology", err)
			ids[i] = ti
			partition[ti] = true
		} else {
			ti, err := tre.TipIndex(tq.set2[i-2])
			tipIndexPanic("testQuartet topology", err)
			ids[i] = ti
		}
	}
	slices.Sort(ids)
	var power uint8 = 1
	var topo uint8 = 0b0000
	for _, id := range ids {
		if partition[id] {
			topo |= power
		}
		power *= 2
	}
	return topo
}

func testQuartetEqual(q *Quartet, tq *TestQuartet, tre *tree.Tree) bool {
	result := q.topology ^ tq.topology(tre)
	return result == 0b0000 || result == 0b1111
}

func (tq *TestQuartet) String() string {
	return fmt.Sprintf("%s%s|%s%s", tq.set1[0], tq.set1[1], tq.set2[0], tq.set2[1])
}

func setToString(qSet map[Quartet]bool, tre *tree.Tree) string {
	str := "{"
	for q := range qSet {
		str += q.String(tre) + ", "
	}
	return str + "}"
}

func listToString(qList []*Quartet, tre *tree.Tree) string {
	str := "{"
	for _, q := range qList {
		str += q.String(tre) + ", "
	}
	return str + "}"
}

func stringListToQMap(list []string, tre *tree.Tree) (map[Quartet]bool, error) {
	qSet := make(map[Quartet]bool)
	for _, nwk := range list {
		tr, err := newick.NewParser(strings.NewReader(nwk)).Parse()
		if err != nil {
			return nil, fmt.Errorf("invalid newick tree %s; test is written wrong", nwk)
		}
		q, err := NewQuartet(tr, tre)
		if err != nil {
			return nil, fmt.Errorf("invalid newick tree %s; test is written wrong", nwk)
		}
		qSet[*q] = true
	}
	return qSet, nil
}
