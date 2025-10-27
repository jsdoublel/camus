package infer

import (
	"os"
	"runtime"
	"strings"
	"testing"

	"github.com/evolbioinfo/gotree/io/newick"
	"github.com/evolbioinfo/gotree/tree"

	gr "github.com/jsdoublel/camus/internal/graphs"
	pr "github.com/jsdoublel/camus/internal/prep"
	sc "github.com/jsdoublel/camus/internal/score"
)

func TestInfer(t *testing.T) {
	testCases := []struct {
		name        string
		constTree   string
		geneTrees   []string
		expNumEdges int
		result      string
	}{
		{
			name:      "basic one-edge",
			constTree: "(A,(B,(C,(D,(E,(F,(G,(H,(I,J)))))))));",
			geneTrees: []string{
				"(A,(B,(C,D)));",
				"(B,(C,D),E);",
			},
			expNumEdges: 1,
			result:      "(A,(B,((C)#H1,((#H1,D),(E,(F,(G,(H,(I,J)))))))));",
		},
		{
			name:      "basic two-edges",
			constTree: "((A,((((B,C),D),E),F)),(G,H));",
			geneTrees: []string{
				"((A,B),(C,D));",
				"((G,F),(A,H));",
			},
			expNumEdges: 2,
			result:      "(((A)#H1,((((B,(C)#H2),(#H2,D)),E),F)),(G,(#H1,H)));",
		},
		{
			name:      "two-edge two",
			constTree: "(A,(B,(C,(D,(E,(F,(G,(H,(I,J)))))))));",
			geneTrees: []string{
				"((J,G),(H,I));",
				"((C,G),(E,F));",
			},
			expNumEdges: 2,
			result:      "(A,(B,(C,(D,((E)#H1,((#H1,F),(G,((H)#H2,((#H2,I),J)))))))));",
		},
		{
			name:      "two-edge case two",
			constTree: "((A,((((B,C),D),E),F)),(G,H));",
			geneTrees: []string{
				"((A,B),(C,D));",
				"((A,F),(G,E));",
			},
			expNumEdges: 2,
			result:      "(((A)#H1,((((B,(C)#H2),(#H2,D)),E),(#H1,F))),(G,H));",
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
			expNumEdges: 1,
			result:      "((#H1,((((A)#H1,B),C),D)),E);",
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
			expNumEdges: 2,
			result:      "((#H1,(((((((#H2,((((A)#H2,B),C),D)))#H1,E),F),G),H),I)),J);",
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
			expNumEdges: 1,
			result:      "(((#H1,((((A)#H1,B),C),D)),E),F);",
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
			expNumEdges: 1,
			result:      "((#H1,(((((A,B))#H1,C),D),E)),F);",
		},
		{
			name:      "avoid over-adding edges",
			constTree: "(R,((A,(((B,C),D),((E,F),G))),H));",
			geneTrees: []string{
				"((C,D),(B,H));",
				"((F,G),(E,H));",
				"((R,A),(B,H));",
			},
			expNumEdges: 2,
			result:      "(R,((A,((((B)#H1,C),D),((E,(F)#H2),(#H2,G)))),(#H1,H)));",
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
			expNumEdges: 2,
			result:      "(R,((A,(((B,(C)#H2),(#H2,D)),(((#H1,E),F),G))),(H)#H1));",
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
			expNumEdges: 3,
			result:      "(R,(((A)#H1,(I,(#H1,J))),(((#H2,((B,(C)#H3),(#H3,D))),H),(((E)#H2,F),G))));",
		},
		{
			name:      "cycle below base of one-sided cycle",
			constTree: "(((((A,B),C),(D,(F,(G,H)))),E),R);",
			geneTrees: []string{
				"((B,C),(D,A));",
				"((B,C),(D,E));",
				"((B,C),(A,E));",
				"((B,D),(A,E));",
				"((C,D),(A,E));",
				"((F,G),(H,E));",
				"((F,G),(H,E));",
				"((E,C),(A,R));",
			},
			expNumEdges: 2,
			result:      "((#H1,(((((A)#H1,B),C),(D,((F)#H2,((#H2,G),H)))),E)),R);",
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
		qopts, _ := pr.SetQuartetFilterOptions(0, 0)
		results, err := Infer(constTree, geneTrees, InferOptions{runtime.GOMAXPROCS(0), qopts, &sc.MaximizeScorer{}, false, 0})
		if err != nil {
			t.Fatalf("Infer failed with error %s", err)
		}
		if len(results.Branches) != test.expNumEdges {
			t.Errorf("inferred number of edges %d not equal to expected %d", len(results.Branches), test.expNumEdges)
		}
		for i, res := range results.Branches {
			if len(res) != i+1 {
				t.Errorf("unexpected number of branches %d, expected %d", len(res), i+1)
			}
		}
		result := gr.MakeNetwork(results.Tree, results.Branches[len(results.Branches)-1]).Newick()
		if result != test.result {
			t.Errorf("result %s != expected %s", result, test.result)
		}
	}
}

func TestInfer_Large(t *testing.T) {
	testCases := []struct {
		name          string
		constTreeFile string
		geneTreesFile string
		qMode         int
		filter        float64
		scorer        sc.InitableScorer
		alpha         float64
		expNumEdges   int
		resultFile    string
	}{
		{
			name:          "pauls data default",
			constTreeFile: "testdata/constraint.nwk",
			geneTreesFile: "testdata/gene-trees.nwk",
			qMode:         0,
			filter:        0,
			scorer:        &sc.MaximizeScorer{},
			alpha:         0,
			expNumEdges:   5,
			resultFile:    "testdata/network.nwk",
		},
		{
			name:          "pauls data quartet filter",
			constTreeFile: "testdata/constraint.nwk",
			geneTreesFile: "testdata/gene-trees.nwk",
			qMode:         2,
			filter:        0.5,
			scorer:        &sc.MaximizeScorer{},
			alpha:         0,
			expNumEdges:   4,
			resultFile:    "testdata/net_q2_t05_max.nwk",
		},
		{
			name:          "pauls data norm",
			constTreeFile: "testdata/constraint.nwk",
			geneTreesFile: "testdata/gene-trees.nwk",
			qMode:         2,
			filter:        0.5,
			scorer:        &sc.NormalizedScorer{},
			alpha:         0,
			expNumEdges:   4,
			resultFile:    "testdata/net_q2_t05_norm.nwk",
		},
		{
			name:          "pauls data sym",
			constTreeFile: "testdata/constraint.nwk",
			geneTreesFile: "testdata/gene-trees.nwk",
			qMode:         2,
			filter:        0.5,
			scorer:        &sc.SymDiffScorer{},
			alpha:         0.1,
			expNumEdges:   4,
			resultFile:    "testdata/net_q2_t05_sym_a01.nwk",
		},
	}
	for _, test := range testCases {
		t.Run(test.name, func(t *testing.T) {
			t.Log(test.name)
			inferOpts := BuildTestInferOpts(t, test.qMode, test.filter, test.scorer, test.alpha)
			tre, quartets, err := pr.ReadInputFiles(test.constTreeFile, test.geneTreesFile, pr.Newick)
			if err != nil {
				t.Fatalf("Could not read input files for benchmark (error %s)", err)
			}
			results, err := Infer(tre, quartets.Trees, inferOpts)
			if err != nil {
				t.Fatalf("failed with unexpected err %s", err)
			}
			resultNwks := make([]string, len(results.Branches))
			for i, branches := range results.Branches {
				resultNwks[i] = gr.MakeNetwork(results.Tree, branches).Newick()
			}
			expNwksStr, err := os.ReadFile(test.resultFile)
			if err != nil {
				t.Fatalf("failed with unexpected err %s", err)
			}
			if len(results.Branches) != test.expNumEdges {
				t.Errorf("inferred number of edges %d not equal to expected %d", len(results.Branches), test.expNumEdges)
			}
			for i, res := range results.Branches {
				if len(res) != i+1 {
					t.Errorf("unexpected number of branches %d, expected %d", len(res), i+1)
				}
			}
			expNwks := strings.Split(string(expNwksStr), "\n")
			for i, expNwk := range expNwks[:len(expNwks)-1] {
				if i < len(resultNwks) && strings.TrimSpace(expNwk) != resultNwks[i] {
					t.Errorf("%s != %s, result != expected", resultNwks[i], expNwk)
				}
			}
		})
	}
}

func BuildTestInferOpts(t *testing.T, qmode int, filter float64, scorer sc.InitableScorer, alpha float64) InferOptions {
	t.Helper()
	qopts, err := pr.SetQuartetFilterOptions(qmode, filter)
	if err != nil {
		t.Fatalf("unexpected error while setting quartet filter options")
	}
	return InferOptions{
		NProcs:      runtime.GOMAXPROCS(0),
		QuartetOpts: qopts,
		ScoreMode:   scorer,
		Alpha:       alpha,
	}
}

func BenchmarkInfer(b *testing.B) {
	constTreeFile := "testdata/constraint.nwk"
	geneTreeFile := "testdata/gene-trees.nwk"
	tre, quartets, err := pr.ReadInputFiles(constTreeFile, geneTreeFile, pr.Newick)
	if err != nil {
		b.Fatalf("Could not read input files for benchmark (error %s)", err)
	}
	for b.Loop() {
		qopts, _ := pr.SetQuartetFilterOptions(0, 0)
		_, err := Infer(tre, quartets.Trees, InferOptions{runtime.GOMAXPROCS(0), qopts, &sc.MaximizeScorer{}, false, 0})
		if err != nil {
			b.Fatalf("Infer failed with error %s", err)
		}
	}
}
