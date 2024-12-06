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
			result: "(a,(b,((c)#0,((#0,d),(e,(f,(g,(h,(i,j)))))))));",
		},
		{
			name:      "basic two-edges",
			constTree: "((A,((((B,C),D),E),F)),(G,H));",
			geneTrees: []string{
				"((A,B),(C,D));",
				"((G,F),(A,H));",
			},
			result: "(((A)#0,((((B,(C)#1),(#1,D)),E),F)),(G,(#0,H)));",
		},
		{
			name:      "two-edge two",
			constTree: "(a,(b,(c,(d,(e,(f,(g,(h,(i,j)))))))));",
			geneTrees: []string{
				"((j,g),(h,i));",
				"((c,g),(e,f));",
			},
			result: "(a,(b,(c,(d,((e)#0,((#0,f),(g,((h)#1,((#1,i),j)))))))));",
		},
		{
			name:      "two-edge case two",
			constTree: "((A,((((B,C),D),E),F)),(G,H));",
			geneTrees: []string{
				"((A,B),(C,D));",
				"((A,F),(G,E));",
			},
			result: "(((A)#0,((((B,(C)#1),(#1,D)),E),(#0,F))),(G,H));",
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
