package main

import (
	"flag"
	"fmt"
	"log"
	"os"

	"camus/alg"
	"camus/netio"
)

type args struct {
	treeFile     string
	geneTreeFile string
}

func parseArgs() args {
	flag.Usage = func() {
		fmt.Fprintln(os.Stderr, "usage: camus [-h] <constraint_tree> <gene_tree>")
		fmt.Fprintln(os.Stderr, "\npositional arguments (required):")
		fmt.Fprintln(os.Stderr, "  <constraint_tree>        constraint newick tree")
		fmt.Fprintln(os.Stderr, "  <gene_tree>              gene tree newick file")
		fmt.Fprintln(os.Stderr, "\nflags:")
		flag.PrintDefaults()
		fmt.Fprintln(os.Stderr, "\nexample:")
		fmt.Fprintln(os.Stderr, "  camus contraint.nwk gene-trees.nwk > out.nwk 2> log.txt")
	}
	help := flag.Bool("h", false, "prints this message")
	flag.Parse()
	if *help {
		flag.Usage()
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
