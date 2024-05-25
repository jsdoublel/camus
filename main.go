package main

import (
	"flag"
	"fmt"
	"os"

	"camus/alg"
	"camus/io"
)

type args struct {
	treeFile     string
	quartetsFile string
}

func parseArgs() args {
	flag.NewFlagSet("CAMUS", flag.ContinueOnError)
	treeFile := flag.String("t", "", "constraint tree")
	quartetsFile := flag.String("q", "", "quartets file")
	flag.Parse()
	if *treeFile == "" || *quartetsFile == "" {
		fmt.Fprintln(os.Stderr, "Error: both -t and -q are required")
		flag.Usage()
		os.Exit(1)
	}
	return args{treeFile: *treeFile, quartetsFile: *quartetsFile}
}

func main() {
	args := parseArgs()
	tre, quartets, err := io.ReadInputFiles(args.treeFile, args.quartetsFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error importing file data:\n%+v\n", err)
		os.Exit(2)
	}
	branches, err := alg.CAMUS(tre, quartets)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Sisyphus was not happy :(\n%+v\n", err)
		os.Exit(3)
	}
	fmt.Println(branches)
}
