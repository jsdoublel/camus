package infer

import (
	"strings"
	"testing"

	"github.com/evolbioinfo/gotree/io/newick"
	"github.com/evolbioinfo/gotree/tree"

	"github.com/jsdoublel/camus/graphs"
	"github.com/jsdoublel/camus/prep"
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
			constTree: "(A,(B,(C,(D,(E,(F,(G,(H,(I,J)))))))));",
			geneTrees: []string{
				"(A,(B,(C,D)));",
				"(B,(C,D),E);",
			},
			result: "(A,(B,((C)#H0,((#H0,D),(E,(F,(G,(H,(I,J)))))))));",
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
			constTree: "(A,(B,(C,(D,(E,(F,(G,(H,(I,J)))))))));",
			geneTrees: []string{
				"((J,G),(H,I));",
				"((C,G),(E,F));",
			},
			result: "(A,(B,(C,(D,((E)#H0,((#H0,F),(G,((H)#H1,((#H1,I),J)))))))));",
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
		{
			name:      "duplicate quartet basic",
			constTree: "(((((A,B),C),D),E),F);",
			geneTrees: []string{
				"((E,D),(F,C));",
				"((C,B),(A,D));",
				"((D,C),(A,E));",
				"((D,C),(A,E));",
				"((D,C),(A,E));",
			},
			result: "(((#H0,((((A)#H0,B),C),D)),E),F);",
		},
		{
			name:      "duplicate quartet basic 2",
			constTree: "(((((A,B),C),D),E),F);",
			geneTrees: []string{
				"((E,C),(F,B));",
				"((E,C),(F,B));",
				"((E,C),(F,B));",
				"((C,B),(A,D));",
				"((D,C),(A,E));",
			},
			result: "((#H0,(((((A,B))#H0,C),D),E)),F);",
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
		result := graphs.MakeNetwork(td, edges).Newick()
		if result != test.result {
			t.Errorf("result %s != expected %s", result, test.result)
		}
	}
}

func BenchmarkCAMUS(b *testing.B) {
	constTreeFile := "../testdata/benchmark/constraint.nwk"
	geneTreeFile := "../testdata/benchmark/gene-trees.nwk"
	tre, quartets, err := prep.ReadInputFiles(constTreeFile, geneTreeFile)
	if err != nil {
		b.Fatalf("Could not read input files for benchmark (error %s)", err)
	}
	for i := 0; i < b.N; i++ {
		CAMUS(tre, quartets)
	}
}
