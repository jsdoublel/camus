package io

import (
	"bufio"
	"fmt"
	"os"
	"strings"

	"github.com/evolbioinfo/gotree/io/newick"
	"github.com/evolbioinfo/gotree/tree"
)

func ReadInputFiles(treeFile, quartetsFile string) (*tree.Tree, []*tree.Tree, error) {
	tre, err := readTreeFile(treeFile)
	if err != nil {
		return nil, nil, err
	}
	quartets, err := readQuartetsFile(quartetsFile)
	if err != nil {
		return nil, nil, err
	}
	return tre, quartets, nil
}

func readTreeFile(treeFile string) (*tree.Tree, error) {
	treStr, err := os.ReadFile(treeFile)
	if err != nil {
		return nil, fmt.Errorf("Error reading tree file:\n%w", err)
	}
	if strings.Count(string(treStr), ";") != 1 {
		return nil, fmt.Errorf("There should only be exactly one newick tree in tree file: %s", treeFile)
	}
	tre, err := newick.NewParser(strings.NewReader(string(treStr))).Parse()
	if err != nil {
		return nil, fmt.Errorf("Error parsing tree newick string:\n%w", err)
	}
	return tre, nil
}

func readQuartetsFile(quartetsFile string) ([]*tree.Tree, error) {
	file, err := os.Open(quartetsFile)
	if err != nil {
		return nil, fmt.Errorf("Error opening %s:\n%w", quartetsFile, err)
	}
	defer file.Close()
	scanner := bufio.NewScanner(file)
	quartetList := []*tree.Tree{}
	for i := 0; scanner.Scan(); i++ {
		line := strings.TrimSpace(scanner.Text())
		if line != "" {
			quartet, err := newick.NewParser(strings.NewReader(line)).Parse()
			if err != nil {
				return nil, fmt.Errorf("Error reading quartet on line %d:\n%w", i, err)
			}
			quartetList = append(quartetList, quartet)
		}
	}
	return quartetList, nil
}
