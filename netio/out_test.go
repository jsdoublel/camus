package netio

import (
	"camus/prep"
	"strings"
	"testing"

	"github.com/evolbioinfo/gotree/io/newick"
	"github.com/evolbioinfo/gotree/tree"
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
			result:    "((A,(B,(C,(#0,F))a)b)c,(D,(E)#0)d)e;",
		},
	}
	for _, test := range testCases {
		t.Run(test.name, func(t *testing.T) {
			constTree, err := newick.NewParser(strings.NewReader(test.constTree)).Parse()
			if err != nil {
				t.Fatalf("%s cannot be parsed as newick. Test case is written incorrectly", test.constTree)
			}
			t.Logf("root %s", constTree.Root().Name())
			td, err := prep.Preprocess(constTree, []*tree.Tree{})
			if err != nil {
				t.Fatalf("%s can not be evaluated (err: %s). Test case is written incorrectly", constTree.Newick(), err)
			}
			edges := make([][2]int, len(test.edges))
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
				edges[i] = [2]int{u[0].Id(), w[0].Id()}
			}
			t.Logf("edges %v", edges)
			result := MakeNetwork(td, edges)
			if result != test.result {
				t.Errorf("%s != %s", result, test.result)
			}
		})
	}
}
