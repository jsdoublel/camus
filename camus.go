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

	gr "github.com/jsdoublel/camus/internal/graphs"
	in "github.com/jsdoublel/camus/internal/infer"
	pr "github.com/jsdoublel/camus/internal/prep"
	sc "github.com/jsdoublel/camus/internal/score"
)

const (
	Version    = "v0.7.2"
	ErrMessage = "camus incountered an error ::"

	Infer Command = iota
	Score

	DefaultFormat    = "newick"
	DefaultScoreMode = "max"
	DefaultQMode     = 0
	DefaultThreshold = 0.5
	DefaultAlpha     = 0.1
)

type Command int

var parseCommand = map[string]Command{
	"infer": Infer,
	"score": Score,
}

type args struct {
	command      Command         // infer or score
	gtFormat     pr.Format       // gene tree file format
	treeFile     string          // constraint or network tree file
	geneTreeFile string          // gene trees
	inferOpts    in.InferOptions // camus options
}

func parseArgs() args {
	flag.Usage = func() {
		fmt.Fprint(os.Stderr,
			"usage: camus [flags]... [command] <tree_file> <gene_tree_file>\n",
			"\n",
			"commands:\n\n",
			"  infer\t\t\tfinds level-1 networks given constraint tree and gene trees\n",
			"  score\t\t\tscore each reticulation branch with respects to gene trees\n",
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
			"  infer command example:\n",
			"\tcamus infer constraint.nwk gene-trees.nwk > network.nwk 2> log.txt\n\n",
			"  score command example:\n",
			"\tcamus score network.nwk gene-trees.nwk > scores.csv 2> log.txt\n",
		)
	}
	format, ok := pr.ParseFormat[DefaultFormat]
	if !ok {
		panic(fmt.Sprintf("bad default format %s", DefaultFormat))
	}
	flag.Var(&format, "f", "gene tree `format` [newick|nexus] (default \"newick\")")
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
	if flag.NArg() != 3 {
		parserError("three positional arguments required: <command> <tree> <gene_tree_file>")
	}
	cmd, ok := parseCommand[flag.Arg(0)]
	if !ok {
		parserError(fmt.Sprintf("\"%s\" is not a valid command: either \"infer\" or \"score\" required", flag.Arg(0)))
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
		command:      cmd,
		gtFormat:     format,
		treeFile:     flag.Arg(1),
		geneTreeFile: flag.Arg(2),
		inferOpts:    *inferOpts,
	}
}

// prints message, usage, and exits (statis code 1)
func parserError(message string) {
	fmt.Fprintln(os.Stderr, message+"\n")
	flag.Usage()
	os.Exit(1)
}

func parsePath(binPath string) string {
	parts := strings.Split(binPath, string(os.PathSeparator))
	return parts[len(parts)-1]
}

func main() {
	log.SetFlags(log.LstdFlags | log.Lmicroseconds)
	args := parseArgs()
	log.Printf("camus %s", Version)
	log.Printf("invoked as: %s", strings.Join(append([]string{parsePath(os.Args[0])}, os.Args[1:]...), " "))
	tre, geneTrees, err := pr.ReadInputFiles(args.treeFile, args.geneTreeFile, args.gtFormat)
	if err != nil {
		log.Fatalf("%s %s\n", ErrMessage, err)
	}
	switch args.command {
	case Infer:
		log.Println("running infer...")
		td, results, err := in.Infer(tre, geneTrees.Trees, args.inferOpts)
		if err != nil {
			log.Fatalf("%s %s\n", ErrMessage, err)
		}
		for _, branches := range results {
			fmt.Println(gr.MakeNetwork(td, branches).Newick())
		}
	case Score:
		if !args.inferOpts.QuartetOpts.QuartetFilterOff() {
			log.Println("WARNING: quartet mode != 0 is not supported for score command at this time. Defaulting to 0.")
		}
		log.Println("running score...")
		network, err := pr.ConvertToNetwork(tre)
		if err != nil {
			log.Fatalf("%s %s\n", ErrMessage, err)
		}
		scores, err := sc.ReticulationScore(network, geneTrees.Trees)
		if err != nil {
			log.Fatalf("%s %s\n", ErrMessage, err)
		}
		if err := pr.WriteRetScoresToCSV(scores, geneTrees.Names); err != nil {
			log.Fatalf("%s %s\n", ErrMessage, err)
		}
	default:
		panic(fmt.Sprintf("invalid command (%d)", args.command))
	}
}
