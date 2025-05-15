package infer

import (
	"os"
	"strings"
	"testing"

	"github.com/evolbioinfo/gotree/io/newick"
	"github.com/evolbioinfo/gotree/tree"

	gr "github.com/jsdoublel/camus/graphs"
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
				"((E,I),(A,J));",
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
		{
			name:      "avoid over-adding edges",
			constTree: "(R,((A,(((B,C),D),((E,F),G))),H));",
			geneTrees: []string{
				"((C,D),(B,H));",
				"((F,G),(E,H));",
				"((R,A),(B,H));",
			},
			result: "(R,((A,((((B)#H0,C),D),((E,(F)#H1),(#H1,G)))),(#H0,H)));",
		},
		{
			name:      "avoid over-adding edges 2",
			constTree: "(R,((A,(((B,C),D),((E,F),G))),H));",
			geneTrees: []string{
				"((C,D),(B,H));",
				"((F,G),(E,H));",
				"((R,A),(B,H));",
				"((R,D),(E,H));",
			},
			result: "(R,((A,(((B,(C)#H1),(#H1,D)),(((#H0,E),F),G))),(H)#H0));",
		},
		{
			name:      "test under node u lookup",
			constTree: "(R,((A,(I,J)),((((B,C),D),H),((E,F),G))));",
			geneTrees: []string{
				"((C,D),(B,A));",
				"((F,G),(E,A));",
				"((H,A),(E,C));",
				"((H,R),(E,C));",
				"((I,R),(J,A));",
			},
			result: "(R,(((A)#H0,(I,(#H0,J))),(((#H1,((B,(C)#H2),(#H2,D))),H),(((E)#H1,F),G))));",
		},
	}
	for _, test := range testCases {
		constTree, err := newick.NewParser(strings.NewReader(test.constTree)).Parse()
		if err != nil {
			t.Fatalf("cannot parse %s as newick tree", test.constTree)
		}
		geneTrees := make([]*tree.Tree, len(test.geneTrees))
		for i, g := range test.geneTrees {
			geneTrees[i], err = newick.NewParser(strings.NewReader(g)).Parse()
			if err != nil {
				t.Fatalf("cannot parse %s as newick tree", g)
			}
		}
		td, results, err := CAMUS(constTree, geneTrees)
		if err != nil {
			t.Fatalf("CAMUS failed with error %s", err)
		}
		result := gr.MakeNetwork(td, results[len(results)-1]).Newick()
		if result != test.result {
			t.Errorf("result %s != expected %s", result, test.result)
		}
	}
}

func TestCAMUS_Large(t *testing.T) {
	testCases := []struct {
		name      string
		constTree string
		geneTrees string
		result    string
	}{
		{
			name:      "pauls data",
			constTree: "../testdata/large/constraint.nwk",
			geneTrees: "../testdata/large/gene-trees.nwk",
			result:    "../testdata/large/network.nwk",
		},
	}
	for _, test := range testCases {
		t.Run(test.name, func(t *testing.T) {
			tre, quartets, err := prep.ReadInputFiles(test.constTree, test.geneTrees, prep.Newick)
			if err != nil {
				t.Fatalf("Could not read input files for benchmark (error %s)", err)
			}
			td, results, err := CAMUS(tre, quartets.Trees)
			if err != nil {
				t.Fatalf("failed with unexpected err %s", err)
			}
			ntw := gr.MakeNetwork(td, results[len(results)-1])
			bts, err := os.ReadFile(test.result)
			if err != nil {
				t.Fatalf("failed with unexpected err %s", err)
			}
			if strings.TrimSpace(string(bts)) != ntw.Newick() {
				t.Errorf("%s != %s, result != expected", ntw.Newick(), string(bts))
			}
		})
	}
}

func BenchmarkCAMUS(b *testing.B) {
	constTreeFile := "../testdata/large/constraint.nwk"
	geneTreeFile := "../testdata/large/gene-trees.nwk"
	tre, quartets, err := prep.ReadInputFiles(constTreeFile, geneTreeFile, prep.Newick)
	if err != nil {
		b.Fatalf("Could not read input files for benchmark (error %s)", err)
	}
	for b.Loop() {
		CAMUS(tre, quartets.Trees)
	}
}
