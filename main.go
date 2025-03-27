/*
CAMUS (Constrained Algorithm Maximizing qUartetS) is a dynamic programming
algorithm for inferring level-1 phylogenetic networks from quartets and a
constraint tree.

usage: camus [-h] <constraint_tree> <gene_tree>

	positional arguments (required):
	  <constraint_tree>        constraint newick tree
	  <gene_tree>              gene tree newick file

	flags:
	 -h	prints this message and exits

	example:
	  camus contraint.nwk gene-trees.nwk > out.nwk 2> log.txt
*/
package main

import (
	"flag"
	"fmt"
	"log"
	"os"

	"github.com/jsdoublel/camus/alg"
	"github.com/jsdoublel/camus/netio"
)

var version = "v0.1.2"

type args struct {
	treeFile     string
	geneTreeFile string
}

func parseArgs() args {
	flag.Usage = func() {
		fmt.Fprintln(os.Stderr,
			"usage: camus [-h| -v] <constraint_tree> <gene_trees>\n",
			"\n",
			"positional arguments (required):\n",
			"  <constraint_tree>        constraint newick tree\n",
			"  <gene_trees>             gene tree newick file\n",
			"\n",
			"flags:",
		)
		flag.PrintDefaults()
		fmt.Fprintln(os.Stderr,
			"\n",
			"example:\n",
			"  camus contraint.nwk gene-trees.nwk > out.nwk 2> log.txt",
		)
	}
	help := flag.Bool("h", false, "prints this message and exits")
	ver := flag.Bool("v", false, "prints version number and exits")

	flag.Parse()
	if *help {
		flag.Usage()
		os.Exit(0)
	}
	if *ver {
		fmt.Printf("CAMUS version %s", version)
		os.Exit(0)
	}
	if flag.NArg() != 2 {
		fmt.Fprintln(os.Stderr, "error: two required positional arguments <constraint_tree> <gene_tree_file>")
		flag.Usage()
		os.Exit(1)
	}
	return args{treeFile: flag.Arg(0), geneTreeFile: flag.Arg(1)}
}

func main() {
	log.SetFlags(log.LstdFlags | log.Lmicroseconds)
	args := parseArgs()
	log.Printf("CAMUS version %s", version)
	constTree, geneTrees, err := netio.ReadInputFiles(args.treeFile, args.geneTreeFile)
	if err != nil {
		log.Fatalf("error importing file data -- %s\n", err)
	}
	td, branches, err := alg.CAMUS(constTree, geneTrees)
	if err != nil {
		log.Fatalf("Sisyphus was not happy :( -- %s\n", err)
	}
	fmt.Println(netio.MakeNetwork(td, branches))
}
