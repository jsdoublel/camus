package prep

import (
	"bufio"
	"encoding/csv"
	"errors"
	"fmt"
	"os"
	"slices"
	"strconv"
	"strings"

	"github.com/jsdoublel/camus/graphs"

	"github.com/evolbioinfo/gotree/io/newick"
	"github.com/evolbioinfo/gotree/io/nexus"
	"github.com/evolbioinfo/gotree/tree"
)

var (
	ErrInvalidTreeFile = errors.New("invalid tree file")
	ErrInvalidFormat   = errors.New("invalid format")
)

type GeneTrees struct {
	Trees []*tree.Tree // gene trees
	Names []string     // gene names
}

// Reads in and validates constraint tree and gene tree input files.
// Returns an error if the newick format is invalid, or the file is invalid for
// some other reason (e.g., more than one constraint tree)
func ReadInputFiles(treeFile, genetreesFile, format string) (*tree.Tree, *GeneTrees, error) {
	tre, err := readTreeFile(treeFile)
	if err != nil {
		return nil, nil, err
	}
	genetrees, err := readGeneTreesFile(genetreesFile, format)
	if err != nil {
		return nil, nil, err
	}
	return tre, genetrees, nil
}

// reads and validates constraint tree file
func readTreeFile(treeFile string) (*tree.Tree, error) {
	treBytes, err := os.ReadFile(treeFile)
	if err != nil {
		return nil, fmt.Errorf("error reading tree file: %w", err)
	}
	treStr := strings.TrimSpace(string(treBytes))
	if strings.Count(treStr, "\n") != 0 || treStr == "" {
		return nil, fmt.Errorf("%w, there should only be exactly one newick tree in tree file %s", ErrInvalidTreeFile, treeFile)
	}
	tre, err := newick.NewParser(strings.NewReader(treStr)).Parse()
	if err != nil {
		return nil, fmt.Errorf("%w, error parsing tree newick string from %s: %s", ErrInvalidFormat, treeFile, err.Error())
	}
	return tre, nil
}

// reads and validates gene tree file
func readGeneTreesFile(genetreesFile, format string) (*GeneTrees, error) {
	file, err := os.Open(genetreesFile)
	if err != nil {
		return nil, fmt.Errorf("error opening %s, %w", genetreesFile, err)
	}
	defer file.Close()
	geneTreeList := make([]*tree.Tree, 0)
	geneTreeNames := make([]string, 0)
	switch format {
	case "newick":
		scanner := bufio.NewScanner(file)
		for i := 0; scanner.Scan(); i++ {
			line := strings.TrimSpace(scanner.Text())
			if line != "" {
				genetree, err := newick.NewParser(strings.NewReader(line)).Parse()
				if err != nil {
					return nil, fmt.Errorf("%w, error reading gene tree on line %d in %s: %s", ErrInvalidFormat, i, genetreesFile, err.Error())
				}
				geneTreeList = append(geneTreeList, genetree)
			}
		}
		if len(geneTreeList) < 1 {
			return nil, fmt.Errorf("%w, empty gene tree file %s", ErrInvalidTreeFile, genetreesFile)
		}
		geneTreeNames = make([]string, 0)
		for i := range len(geneTreeList) {
			geneTreeNames = append(geneTreeNames, strconv.Itoa(i+1))
		}
	case "nexus":
		nex, err := nexus.NewParser(file).Parse()
		if err != nil {
			return nil, fmt.Errorf("%w, error reading gene tree nexus file %s: %s", ErrInvalidFormat, genetreesFile, err.Error())
		}
		nex.IterateTrees(func(s string, t *tree.Tree) {
			geneTreeList = append(geneTreeList, t)
			geneTreeNames = append(geneTreeNames, s)
		})
	default:
		return nil, fmt.Errorf("%w, not a valid file format", ErrInvalidTreeFile)
	}
	return &GeneTrees{Trees: geneTreeList, Names: geneTreeNames}, nil
}

// Read in extended newick file and make network
func ConvertToNetwork(ntw *tree.Tree) (network *graphs.Network, err error) {
	if !ntw.Rooted() {
		return nil, fmt.Errorf("%w, network is not rooted", ErrInvalidTree)
	}
	if !NetworkIsBinary(ntw) {
		return nil, fmt.Errorf("%w, network must be fully resolved", ErrInvalidTree)
	}
	ret := make(map[string][2]int)
	defer func() {
		if r := recover(); r != nil {
			err = fmt.Errorf("%w, too many or invalid matching reticulation label %v", ErrInvalidFormat, r)
		}
	}()
	ntw.PostOrder(func(cur, prev *tree.Node, e *tree.Edge) (keep bool) {
		if strings.Contains(cur.Name(), "#") {
			branch := ret[cur.Name()]
			var v *tree.Node
			if cur.Tip() {
				p, err := prev.Parent()
				if err != nil && err.Error() != "The node has no parent : May be the root?" {
					panic(fmt.Sprintf("%s", cur.Name()))
				}
				for _, n := range prev.Neigh() {
					if n != cur && n != p {
						v = n
					}
				}
				if branch[graphs.Ui] != 0 || v == nil {
					panic(fmt.Sprintf("%s", cur.Name()))
				}
				branch[graphs.Ui] = v.Id()
			} else {
				for _, n := range cur.Neigh() {
					if n != cur && n != prev {
						v = n
					}
				}
				if branch[graphs.Wi] != 0 || v == nil {
					panic(fmt.Sprintf("%s", cur.Name()))
				}
				branch[graphs.Wi] = v.Id()
			}
			ret[cur.Name()] = branch
		}
		return true
	})
	if len(ret) == 0 {
		return nil, fmt.Errorf("%w, no reticulations - not a network", ErrInvalidTree)
	}
	for label, branch := range ret {
		if branch[graphs.Ui] == 0 || branch[graphs.Wi] == 0 { // assumes root node is not labeled as reticulation
			return nil, fmt.Errorf("%w, label %s is unmatched", ErrInvalidFormat, label)
		}
	}
	ntw.UpdateTipIndex()
	return &graphs.Network{NetTree: ntw, Reticulations: ret}, nil
}

// Write csv file containing reticulation branch scores to stdout
func WriteBranchScoresToCSV(scores []*map[string]float64, names []string) {
	branchNames := make([]string, 0)
	for k := range *scores[0] {
		branchNames = append(branchNames, k)
	}
	slices.SortFunc(branchNames, func(a, b string) int {
		if diff := len(a) - len(b); diff != 0 {
			return diff
		}
		return strings.Compare(a, b)
	})
	data := make([][]string, len(scores)+1)
	data[0] = append([]string{"gene"}, branchNames...)
	for i, row := range scores {
		data[i+1] = []string{names[i]}
		for _, br := range branchNames {
			data[i+1] = append(data[i+1], strconv.FormatFloat((*row)[br], 'f', -1, 64))
		}
	}
	writer := csv.NewWriter(os.Stdout)
	defer writer.Flush()
	writer.WriteAll(data)
}
