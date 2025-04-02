/*
CAMUS (Constrained Algorithm Maximizing qUartetS) is a dynamic programming
algorithm for inferring level-1 phylogenetic networks from quartets and a
constraint tree.

usage: camus [-h| -v] <constraint_tree> <gene_trees>

positional arguments (required):

	<constraint_tree>        constraint newick tree
	<gene_trees>             gene tree newick file

flags:

	-h	prints this message and exits
	-v	prints version number and exits

example:

	camus contraint.nwk gene-trees.nwk > out.nwk 2> log.txt
*/
package main

import (
	"flag"
	"fmt"
	"log"
	"os"

	"github.com/mattn/go-isatty"

	"github.com/jsdoublel/camus/graphs"
	"github.com/jsdoublel/camus/infer"
	"github.com/jsdoublel/camus/prep"
)

var version = "v0.1.3"

type args struct {
	treeFile     string
	geneTreeFile string
}

var ErrMessage = fmt.Sprintf("[%s]:", red("error"))

func parseArgs() args {
	flag.Usage = func() {
		fmt.Fprint(os.Stderr,
			"usage: camus [-h| -v] <constraint_tree> <gene_trees>\n",
			"\n",
			"positional arguments (required):\n\n",
			"  <constraint_tree>\tconstraint newick tree\n",
			"  <gene_trees>\t\tgene tree newick file\n",
			"\n",
			"flags:\n\n",
		)
		flag.PrintDefaults()
		fmt.Fprint(os.Stderr,
			"\n",
			"example:\n\n",
			"  camus contraint.nwk gene-trees.nwk > out.nwk 2> log.txt\n",
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
		fmt.Fprintf(os.Stderr, "%s two positional arguments required: <constraint_tree> <gene_tree_file>\n", ErrMessage)
		flag.Usage()
		os.Exit(1)
	}
	return args{treeFile: flag.Arg(0), geneTreeFile: flag.Arg(1)}
}

// makes text printed to terminal red
func red(str string) string {
	if isatty.IsTerminal(os.Stderr.Fd()) {
		return fmt.Sprintf("\033[31m%s\033[0m", str)
	}
	return str
}

func main() {
	log.SetFlags(log.LstdFlags | log.Lmicroseconds)
	args := parseArgs()
	log.Printf("CAMUS version %s", version)
	constTree, geneTrees, err := prep.ReadInputFiles(args.treeFile, args.geneTreeFile)
	if err != nil {
		log.Fatalf("%s %s\n", ErrMessage, err)
	}
	td, branches, err := infer.CAMUS(constTree, geneTrees)
	if err != nil {
		log.Fatalf("%s %s\n", ErrMessage, err)
	}
	fmt.Println(graphs.MakeNetwork(td, branches).Newick())
}
