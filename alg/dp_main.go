package alg

import (
	"fmt"

	"github.com/evolbioinfo/gotree/tree"
)

func CAMUS(tre *tree.Tree, rawQuartets []*tree.Tree) ([][]uint, error) {
	quartets, lcaMat, leafsets, err := Preprocess(tre, rawQuartets)
	if err != nil {
		return nil, fmt.Errorf("Preprocess error: %w", err)
	}
	fmt.Println(quartets)
	fmt.Println(lcaMat)
	fmt.Println(leafsets)
	return nil, nil
}
