package alg

import (
	"camus/netio"
	"strings"
	"testing"

	"github.com/evolbioinfo/gotree/io/newick"
	"github.com/evolbioinfo/gotree/tree"
)

func TestCAMUS(t *testing.T) {
	testCases := []struct {
		name      string
		constTree string
		geneTrees []string
		result    string
	}{
		{
			name:      "basic one-edge",
			constTree: "(a,(b,(c,(d,(e,(f,(g,(h,(i,j)))))))));",
			geneTrees: []string{
				"(a,(b,(c,d)));",
				"(b,(c,d),e);",
			},
			result: "(a,(b,((c)#H0,((#H0,d),(e,(f,(g,(h,(i,j)))))))));",
		},
		{
			name:      "basic two-edges",
			constTree: "((A,((((B,C),D),E),F)),(G,H));",
			geneTrees: []string{
				"((A,B),(C,D));",
				"((G,F),(A,H));",
			},
			result: "(((A)#H0,((((B,(C)#H1),(#H1,D)),E),F)),(G,(#H0,H)));",
		},
		{
			name:      "two-edge two",
			constTree: "(a,(b,(c,(d,(e,(f,(g,(h,(i,j)))))))));",
			geneTrees: []string{
				"((j,g),(h,i));",
				"((c,g),(e,f));",
			},
			result: "(a,(b,(c,(d,((e)#H0,((#H0,f),(g,((h)#H1,((#H1,i),j)))))))));",
		},
		{
			name:      "two-edge case two",
			constTree: "((A,((((B,C),D),E),F)),(G,H));",
			geneTrees: []string{
				"((A,B),(C,D));",
				"((A,F),(G,E));",
			},
			result: "(((A)#H0,((((B,(C)#H1),(#H1,D)),E),(#H0,F))),(G,H));",
		},
		{
			name:      "one-sided cycle test",
			constTree: "((((A,B),C),D),E);",
			geneTrees: []string{
				"((B,C),(D,A));",
				"((B,C),(D,E));",
				"((B,C),(A,E));",
				"((B,D),(A,E));",
				"((C,D),(A,E));",
			},
			result: "((#H0,((((A)#H0,B),C),D)),E);",
		},
		{
			name:      "double one-sided cycle test",
			constTree: "(((((((((A,B),C),D),E),F),G),H),I),J);",
			geneTrees: []string{
				"((B,C),(D,A));",
				"((B,C),(D,E));",
				"((B,C),(A,E));",
				"((B,D),(A,E));",
				"((C,D),(A,E));",
				"((G,H),(I,D));",
				"((G,H),(I,J));",
				"((G,H),(D,J));",
				"((G,I),(D,J));",
				"((H,I),(D,J));",
			},
			result: "((#H0,(((((((#H1,((((A)#H1,B),C),D)))#H0,E),F),G),H),I)),J);",
		},
	}
	for _, test := range testCases {
		constTree, err := newick.NewParser(strings.NewReader(test.constTree)).Parse()
		if err != nil {
			t.Errorf("cannot parse %s as newick tree", test.constTree)
		}
		geneTrees := make([]*tree.Tree, len(test.geneTrees))
		for i, g := range test.geneTrees {
			geneTrees[i], err = newick.NewParser(strings.NewReader(g)).Parse()
			if err != nil {
				t.Errorf("cannot parse %s as newick tree", g)
			}
		}
		td, edges, err := CAMUS(constTree, geneTrees)
		if err != nil {
			t.Errorf("CAMUS failed with error %s", err)
		}
		result := netio.MakeNetwork(td, edges)
		if result != test.result {
			t.Errorf("result %s != expected %s", result, test.result)
		}
	}
}
