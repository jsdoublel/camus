package score

import (
    "strings"
    "testing"

    "github.com/evolbioinfo/gotree/io/newick"

    gr "github.com/jsdoublel/camus/internal/graphs"
)

// helper to fetch a node id by label
func nodeIDByLabel(t *testing.T, td *gr.TreeData, label string) int {
    t.Helper()
    nodes, err := td.SelectNodes(label)
    if err != nil {
        t.Fatalf("failed to select node %q: %v", label, err)
    }
    if len(nodes) != 1 {
        t.Fatalf("expected exactly one node for label %q, got %d", label, len(nodes))
    }
    return nodes[0].Id()
}

func TestGetNumTaxaUnderNodes(t *testing.T) {
    testCases := []struct {
        name     string
        newick   string
        uLabel   string
        wLabel   string
        expected []uint
    }{
        {
            name:   "chain_u_b_w_a",
            newick: "((((A,B)a,C)b,D)c,E)r;",
            uLabel: "b",
            wLabel: "a",
            expected: []uint{2, 1, 1, 2}, // [a, sib@b(C), sib@c(D), outside]
        },
        {
            name:   "chain_u_c_w_a",
            newick: "((((A,B)a,C)b,D)c,E)r;",
            uLabel: "c",
            wLabel: "a",
            expected: []uint{2, 1, 1, 1, 1}, // [a, sib@b(C), sib@c(D), sib@r(E), outside]
        },
    }

    for _, tc := range testCases {
        t.Run(tc.name, func(t *testing.T) {
            tre, err := newick.NewParser(strings.NewReader(tc.newick)).Parse()
            if err != nil {
                t.Fatalf("invalid newick in test: %v", err)
            }
            if err := tre.UpdateTipIndex(); err != nil {
                t.Fatalf("failed to update tip index: %v", err)
            }
            td := gr.MakeTreeData(tre, nil)

            u := nodeIDByLabel(t, td, tc.uLabel)
            w := nodeIDByLabel(t, td, tc.wLabel)

            got := getNumTaxaUnderNodes(u, w, td)
            if len(got) != len(tc.expected) {
                t.Fatalf("wrong subsets length: got %d want %d (%v vs %v)", len(got), len(tc.expected), got, tc.expected)
            }
            for i := range got {
                if got[i] != tc.expected[i] {
                    t.Fatalf("subset[%d] = %d, want %d (full=%v)", i, got[i], tc.expected[i], got)
                }
            }
        })
    }
}

func TestCalculatePenalty(t *testing.T) {
    testCases := []struct {
        name     string
        newick   string
        uLabel   string
        wLabel   string
        expected uint
    }{
        {
            name:     "b_to_a",
            newick:   "((((A,B)a,C)b,D)c,E)r;",
            uLabel:   "b",
            wLabel:   "a",
            expected: 4, // subsets [2,1,1,2] -> 2 * (1*1*2)
        },
        {
            name:     "c_to_a",
            newick:   "((((A,B)a,C)b,D)c,E)r;",
            uLabel:   "c",
            wLabel:   "a",
            expected: 8, // subsets [2,1,1,1,1] -> 2 * (4 choose triples of 1s = 4)
        },
    }

    for _, tc := range testCases {
        t.Run(tc.name, func(t *testing.T) {
            tre, err := newick.NewParser(strings.NewReader(tc.newick)).Parse()
            if err != nil {
                t.Fatalf("invalid newick in test: %v", err)
            }
            if err := tre.UpdateTipIndex(); err != nil {
                t.Fatalf("failed to update tip index: %v", err)
            }
            td := gr.MakeTreeData(tre, nil)

            u := nodeIDByLabel(t, td, tc.uLabel)
            w := nodeIDByLabel(t, td, tc.wLabel)

            got := calculatePenalty(u, w, td)
            if got != tc.expected {
                t.Fatalf("penalty = %d, want %d", got, tc.expected)
            }
        })
    }
}
