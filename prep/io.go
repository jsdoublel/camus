package prep

import (
	"bufio"
	"encoding/csv"
	"errors"
	"fmt"
	"os"
	"strconv"
	"strings"

	"github.com/jsdoublel/camus/graphs"

	"github.com/evolbioinfo/gotree/io/newick"
	"github.com/evolbioinfo/gotree/tree"
)

var (
	ErrInvalidTreeFile = errors.New("invalid tree file")
	ErrInvalidNewick   = errors.New("invalid newick format")
)

// Reads in and validates constraint tree and gene tree input files.
// Returns an error if the newick format is invalid, or the file is invalid for
// some other reason (e.g., more than one constraint tree)
func ReadInputFiles(treeFile, genetreesFile string) (*tree.Tree, []*tree.Tree, error) {
	tre, err := readTreeFile(treeFile)
	if err != nil {
		return nil, nil, err
	}
	genetrees, err := readGeneTreesFile(genetreesFile)
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
		return nil, fmt.Errorf("%w, error parsing tree newick string from %s: %s", ErrInvalidNewick, treeFile, err.Error())
	}
	return tre, nil
}

// reads and validates gene tree file
func readGeneTreesFile(genetreesFile string) ([]*tree.Tree, error) {
	file, err := os.Open(genetreesFile)
	if err != nil {
		return nil, fmt.Errorf("error opening %s, %w", genetreesFile, err)
	}
	defer file.Close()
	scanner := bufio.NewScanner(file)
	genetreeList := []*tree.Tree{}
	for i := 0; scanner.Scan(); i++ {
		line := strings.TrimSpace(scanner.Text())
		if line != "" {
			genetree, err := newick.NewParser(strings.NewReader(line)).Parse()
			if err != nil {
				return nil, fmt.Errorf("%w, error reading gene tree on line %d in %s: %s", ErrInvalidNewick, i, genetreesFile, err.Error())
			}
			genetreeList = append(genetreeList, genetree)
		}
	}
	if len(genetreeList) < 1 {
		return nil, fmt.Errorf("%w, empty gene tree file %s", ErrInvalidTreeFile, genetreesFile)
	}
	return genetreeList, nil
}

// Read in extended newick file and make network
func ConvertToNetwork(ntw *tree.Tree) (network *graphs.Network, err error) {
	ret := make(map[string][2]int)
	defer func() {
		if r := recover(); r != nil {
			err = fmt.Errorf("%w, too many or invalid matching reticulation label %v", ErrInvalidNewick, r)
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
	for label, branch := range ret {
		if branch[graphs.Ui] == 0 || branch[graphs.Wi] == 0 { // assumes root node is not labeled as reticulation
			return nil, fmt.Errorf("%w, label %s is unmatched", ErrInvalidNewick, label)
		}
	}
	ntw.UpdateTipIndex()
	return &graphs.Network{NetTree: ntw, Reticulations: ret}, nil
}

// Write csv file containing reticulation branch scores to stdout
func WriteBranchScoresToCSV(scores []*map[string]float64) error {
	header := []string{"gene tree"}
	data := make([][]string, len(scores))
	for k := range *scores[0] {
		header = append(header, k)
	}
	for i, row := range scores {
		data[i] = []string{strconv.Itoa(i)}
		for _, v := range *row {
			data[i] = append(data[i], strconv.FormatFloat(v, 'f', -1, 64))
		}
	}
	writer := csv.NewWriter(os.Stdout)
	defer writer.Flush()
	writer.Write(header)
	writer.WriteAll(data)
	return nil
}
