package main

import (
	"flag"
	"fmt"
	"os"

	"camus/alg"
	"camus/netio"
)

type args struct {
	treeFile     string
	geneTreeFile string
}

func parseArgs() args {
	flag.NewFlagSet("CAMUS", flag.ContinueOnError)
	treeFile := flag.String("t", "", "constraint tree")
	geneTreeFile := flag.String("g", "", "gene tree file")
	flag.Parse()
	if *treeFile == "" || *geneTreeFile == "" {
		fmt.Fprintln(os.Stderr, "Error: both -t and -g are required")
		flag.Usage()
		os.Exit(1)
	}
	return args{treeFile: *treeFile, geneTreeFile: *geneTreeFile}
}

func main() {
	args := parseArgs()
	tre, quartets, err := netio.ReadInputFiles(args.treeFile, args.geneTreeFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error importing file data:\n%+v\n", err)
		os.Exit(2)
	}
	td, branches, err := alg.CAMUS(tre, quartets)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Sisyphus was not happy :(\n%+v\n", err)
		os.Exit(3)
	}
	fmt.Fprintf(os.Stdout, netio.MakeNetwork(td, branches))
}
