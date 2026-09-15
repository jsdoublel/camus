package graphs

import (
	"strings"
	"testing"

	"github.com/evolbioinfo/gotree/io/newick"
)

func TestMakeNetwork(t *testing.T) {
	testCases := []struct {
		name      string
		constTree string
		edges     [][2]string
		result    string
	}{
		{
			name:      "basic",
			constTree: "[&R]((A,(B,(C,F)a)b)c,(D,E)d)e;",
			edges:     [][2]string{{"F", "E"}},
			result:    "((A,(B,(C,(#H1,F))a)b)c,(D,(E)#H1)d)e;",
		},
		{
			name:      "transitive comparator",
			constTree: "[&R]((((A,B)u1,C)mid,G)u3,((D,E)u2,H)mid2)r;",
			edges:     [][2]string{{"u1", "A"}, {"u2", "D"}, {"u3", "u1"}},
			result:    "((#H1,((((#H3,((A)#H3,B)u1))#H1,C)mid,G)u3),((#H2,((D)#H2,E)u2),H)mid2)r;",
		},
	}
	for _, test := range testCases {
		t.Run(test.name, func(t *testing.T) {
			constTree, err := newick.NewParser(strings.NewReader(test.constTree)).Parse()
			if err != nil {
				t.Fatalf("%s cannot be parsed as newick. Test case is written incorrectly", test.constTree)
			}
			t.Logf("root %s", constTree.Root().Name())
			err = constTree.UpdateTipIndex()
			if err != nil {
				t.Error(err)
			}
			td := MakeTreeData(constTree, nil)
			edges := make([]Branch, len(test.edges))
			for i, edge := range test.edges {
				u, err := constTree.SelectNodes(edge[0])
				t.Logf("u: %s, u (id): %d", u[0].Name(), u[0].Id())
				if err != nil || len(u) != 1 {
					t.Fatalf("cannot find node %s or found too many", edge[0])
				}
				w, err := constTree.SelectNodes(edge[1])
				t.Logf("w: %s, w (id): %d", w[0].Name(), w[0].Id())
				if err != nil || len(w) != 1 {
					t.Fatalf("cannot find node %s or found too many", edge[1])
				}
				edges[i] = Branch{IDs: [2]int{u[0].Id(), w[0].Id()}}
			}
			t.Logf("edges %v", edges)
			result := MakeNetwork(td, edges).Newick()
			if result != test.result {
				t.Errorf("%s != %s", result, test.result)
			}
		})
	}
}
