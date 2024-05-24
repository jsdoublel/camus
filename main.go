package main

import (
	"flag"
	"fmt"
	"os"

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
	fmt.Println(tre.AllTipNames())
	fmt.Println(len(quartets))
	// for _, q := range quartets {
	// 	fmt.Printf("%d%d|%d%d\n", q.T1, q.T2, q.T3, q.T4)
	// }
}
