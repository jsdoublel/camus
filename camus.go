/*
CAMUS (Constrained Algorithm Maximizing qUartetS) is a dynamic programming
algorithm for inferring level-1 phylogenetic networks from quartets and a
constraint tree.

usage: camus [flags]... [command] <tree_file> <gene_tree_file>

commands:

	infer			finds level-1 networks given constraint tree and gene trees
	score			score each reticulation branch with respects to gene trees

positional arguments:

	<tree_file>		constraint newick tree (infer) or network (score)
	<gene_tree_file>	gene tree newick file

flags:

	-a float
	  	parameter to adjust penalty for "sym" score mode (default 0.1)
	-asSet
	  	quartet count is calculated as a set (one point per unique topology)
	-f format
	  	gene tree format [newick|nexus] (default "newick")
	-h	prints this message and exits
	-n int
	  	number of parallel processes
	-q int
	  	quartet filter mode number [0, 2] (default 0)
	-s mode
	  	score mode [max|norm|sym] (default "max")
	-t float
	  	threshold for quartet filter [0, 1] (default 0.5)
	-v	prints version number and exits

examples:

	  infer command example:
		camus infer constraint.nwk gene-trees.nwk > network.nwk 2> log.txt

	  score command example:
		camus score network.nwk gene-trees.nwk > scores.csv 2> log.txt
*/
package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"strings"
	"time"

	gr "github.com/jsdoublel/camus/internal/graphs"
	in "github.com/jsdoublel/camus/internal/infer"
	pr "github.com/jsdoublel/camus/internal/prep"
	sc "github.com/jsdoublel/camus/internal/score"
)

const (
	Version    = "v0.8.0"
	ErrMessage = "camus incountered an error ::"
	TimeFormat = "2006-01-02_15-04-05"

	DefaultFormat    = "newick"
	DefaultScoreMode = "max"
	DefaultQMode     = 0
	DefaultThreshold = 0.5
	DefaultAlpha     = 0.1
)

type args struct {
	prefix       string          // outptu prefix
	gtFormat     pr.Format       // gene tree file format
	treeFile     string          // constraint or network tree file
	geneTreeFile string          // gene trees
	inferOpts    in.InferOptions // camus options
}

func parseArgs() args {
	flag.Usage = func() {
		fmt.Fprint(os.Stderr,
			"usage: camus [flags]... <tree_file> <gene_tree_file>\n",
			"\n",
			"positional arguments:\n\n",
			"  <tree_file>\t\tconstraint newick tree (infer) or network (score)\n",
			"  <gene_tree_file>\tgene tree newick file\n",
			"\n",
			"flags:\n\n",
		)
		flag.PrintDefaults()
		fmt.Fprint(os.Stderr,
			"\n",
			"examples:\n\n",
			"\tcamus constraint.nwk gene-trees.nwk > network.nwk 2> log.txt\n\n",
		)
	}
	format, ok := pr.ParseFormat[DefaultFormat]
	if !ok {
		panic(fmt.Sprintf("bad default format %s", DefaultFormat))
	}
	flag.Var(&format, "f", "gene tree `format` [newick|nexus] (default \"newick\")")
	prefix := flag.String("o", "", "output prefix")
	scoreMode := flag.String("s", DefaultScoreMode, "score `mode` [max|norm|sym]")
	mode := flag.Int("q", DefaultQMode, "quartet filter mode number [0, 2] (default 0)")
	thresh := flag.Float64("t", DefaultThreshold, "threshold for quartet filter [0, 1]")
	alpha := flag.Float64("a", DefaultAlpha, "parameter to adjust penalty for \"sym\" score mode, from (0, 1]")
	asSet := flag.Bool("asSet", false, "quartet count is calculated as a set (one point per unique topology)")
	help := flag.Bool("h", false, "prints this message and exits")
	ver := flag.Bool("v", false, "prints version number and exits")
	nprocs := flag.Int("n", 0, "number of parallel processes")
	flag.Parse()
	if *help {
		flag.Usage()
		os.Exit(0)
	}
	if *ver {
		fmt.Printf("camus %s\n", Version)
		os.Exit(0)
	}
	if flag.NArg() != 2 {
		parserError("three positional arguments required: <tree> <gene_tree_file>")
	}
	scorer, ok := sc.ParseScorer[*scoreMode]
	if !ok {
		parserError(fmt.Sprintf("\"%s\" is not a valid score mode: valid score modes are \"max\", \"norm\", and \"sym\"", *scoreMode))
	}
	qOpts, err := pr.SetQuartetFilterOptions(*mode, *thresh)
	if err != nil {
		parserError(err.Error())
	}
	inferOpts, err := in.MakeInferOptions(*nprocs, qOpts, scorer, *asSet, *alpha)
	if err != nil {
		parserError(err.Error())
	}
	return args{
		prefix:       *prefix,
		gtFormat:     format,
		treeFile:     flag.Arg(0),
		geneTreeFile: flag.Arg(1),
		inferOpts:    *inferOpts,
	}
}

// prints message, usage, and exits (status code 1)
func parserError(message string) {
	fmt.Fprintln(os.Stderr, message+"\n")
	flag.Usage()
	os.Exit(1)
}

func defaultPrefix() string {
	parseName := func(s string) string {
		parts := strings.Split(s, string(os.PathSeparator))
		parts = strings.Split(parts[len(parts)-1], ".")
		return strings.Join(parts[:len(parts)-1], ".")
	}
	inputs := fmt.Sprintf("%s_%s", parseName(flag.Arg(0)), parseName(flag.Arg(1)))
	return fmt.Sprintf("camus_%s_%s", inputs, time.Now().Local().Format(TimeFormat))
}

func main() {
	if err := run(); err != nil {
		log.Fatalf("%s %s", ErrMessage, err)
	}
}

func run() error {
	log.SetFlags(log.LstdFlags | log.Lmicroseconds)
	args := parseArgs()
	log.Printf("camus %s", Version)
	log.Printf("invoked as: camus %s", strings.Join(os.Args[1:], " "))
	if args.prefix == "" {
		args.prefix = defaultPrefix()
		log.Printf("output prefix was not set, using \"%s\"", args.prefix)
	}
	tre, geneTrees, err := pr.ReadInputFiles(args.treeFile, args.geneTreeFile, args.gtFormat)
	if err != nil {
		return err
	}
	results, err := in.Infer(tre, geneTrees.Trees, args.inferOpts)
	if err != nil {
		return err
	}
	newicks := make([]string, len(results.Branches))
	for i, branches := range results.Branches {
		newicks[i] = gr.MakeNetwork(results.Tree, branches).Newick()
	}
	if err = pr.WriteDPResultsToCSV(results.Tree, newicks, results.QSatScore, os.Stdout); err != nil {
		return err
	}
	f, err := os.Create(fmt.Sprintf("%s.csv", args.prefix))
	if err != nil {
		return err
	}
	defer func() {
		closeErr := f.Close()
		if closeErr != nil {
			log.Printf("error closing %s.csv, %s", args.prefix, closeErr)
		}
	}()
	if err = pr.WriteDPResultsToCSV(results.Tree, newicks, results.QSatScore, f); err != nil {
		return err
	}
	if err = pr.WriteResultsLineplot(results.QSatScore, args.prefix); err != nil {
		return err
	}
	return nil
}
