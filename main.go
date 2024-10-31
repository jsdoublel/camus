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
	outputFile   string
}

func parseArgs() args {
	flag.NewFlagSet("CAMUS", flag.ContinueOnError)
	treeFile := flag.String("t", "", "constraint tree")
	geneTreeFile := flag.String("g", "", "gene tree file")
	outputFile := flag.String("o", "", "output extended newick file")
	flag.Parse()
	if *treeFile == "" || *geneTreeFile == "" || *outputFile == "" {
		fmt.Fprintln(os.Stderr, "Error: both -t, -g, and -o are required")
		flag.Usage()
		os.Exit(1)
	}
	return args{treeFile: *treeFile, geneTreeFile: *geneTreeFile, outputFile: *outputFile}
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
	out := []byte(netio.MakeNetwork(td, branches))
	err = os.WriteFile(args.outputFile, out, 0644)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error writing output:\n%+v\n", err)
		os.Exit(4)
	}
}
