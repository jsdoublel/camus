package prep

import (
	"bufio"
	"errors"
	"fmt"
	"os"
	"strings"

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
		return nil, fmt.Errorf("error reading constraint tree file: %w", err)
	}
	treStr := strings.TrimSpace(string(treBytes))
	if strings.Count(treStr, "\n") != 0 || treStr == "" {
		return nil, fmt.Errorf("%w, there should only be exactly one newick tree in constraint tree file %s", ErrInvalidTreeFile, treeFile)
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
