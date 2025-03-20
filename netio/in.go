package netio

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
	ErrInvalidConstTree = errors.New("invalid constraint tree")
	ErrInvalidNewick    = errors.New("invalid newick format")
)

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

func readTreeFile(treeFile string) (*tree.Tree, error) {
	treStr, err := os.ReadFile(treeFile)
	if err != nil {
		return nil, fmt.Errorf("error reading constraint tree file: %w", err)
	}
	if strings.Count(string(treStr), ";") != 1 {
		return nil, fmt.Errorf("%w, there should only be exactly one newick tree in constraint tree file %s", ErrInvalidConstTree, treeFile)
	}
	tre, err := newick.NewParser(strings.NewReader(string(treStr))).Parse()
	if err != nil {
		return nil, fmt.Errorf("%w, error parsing tree newick string from %s: %s", ErrInvalidNewick, treeFile, err.Error())
	}
	return tre, nil
}

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
	return genetreeList, nil
}
