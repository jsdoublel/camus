/*
CAMUS (Constrained Algorithm Maximizing qUartetS) is a dynamic programming
algorithm for inferring level-1 phylogenetic networks from quartets and a
constraint tree.

usage: camus [ -f <format> | -h | -v ] <command> <tree> <gene_trees>

commands:

	infer		find level-1 network given constraint tree and gene trees
	score		score each reticulation branch with respects to gene trees

positional arguments:

	<tree>			constraint newick tree (infer) or network (score)
	<gene_trees>	gene tree newick file

flags:

	-f format
	  	gene tree format [ newick | nexus ] (default "newick")
	-h	prints this message and exits
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

	gr "github.com/jsdoublel/camus/graphs"
	"github.com/jsdoublel/camus/infer"
	pr "github.com/jsdoublel/camus/prep"
	"github.com/jsdoublel/camus/score"
)

const (
	Version    = "v0.3.0"
	ErrMessage = "Sisyphus was not happy :("

	Infer Command = iota
	Score
)

type Command int

var parseCommand = map[string]Command{
	"infer": Infer,
	"score": Score,
}

type args struct {
	command      Command   // infer or score
	gtFormat     pr.Format // gene tree file format
	treeFile     string    // constraint or network tree file
	geneTreeFile string    // gene trees
}

func parseArgs() args {
	flag.Usage = func() {
		fmt.Fprint(os.Stderr,
			"usage: camus [ -f <format> | -h | -v ] <command> <tree> <gene_trees>\n",
			"\n",
			"commands:\n\n",
			"  infer\t\tfind level-1 network given constraint tree and gene trees\n",
			"  score\t\tscore each reticulation branch with respects to gene trees\n",
			"\n",
			"positional arguments:\n\n",
			"  <tree>\tconstraint newick tree (infer) or network (score)\n",
			"  <gene_trees>\tgene tree newick file\n",
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
	format := pr.Newick
	flag.Var(&format, "f", "gene tree `format` [ newick | nexus ]")
	help := flag.Bool("h", false, "prints this message and exits")
	ver := flag.Bool("v", false, "prints version number and exits")
	flag.Parse()
	if *help {
		flag.Usage()
		os.Exit(0)
	}
	if *ver {
		fmt.Printf("CAMUS version %s\n", Version)
		os.Exit(0)
	}
	if flag.NArg() != 3 {
		parserError("three positional arguments required: <command> <tree> <gene_tree_file>")
	}
	cmd, ok := parseCommand[flag.Arg(0)]
	if !ok {
		parserError(fmt.Sprintf("\"%s\" is not a valid command: either \"infer\" or \"score\" required", flag.Arg(0)))
	}
	return args{
		command:      cmd,
		gtFormat:     format,
		treeFile:     flag.Arg(1),
		geneTreeFile: flag.Arg(2),
	}
}

// prints message, usage, and exits (statis code 1)
func parserError(message string) {
	fmt.Fprintln(os.Stderr, message)
	flag.Usage()
	os.Exit(1)
}

func main() {
	log.SetFlags(log.LstdFlags | log.Lmicroseconds)
	args := parseArgs()
	log.Printf("CAMUS version %s", Version)
	tre, geneTrees, err := pr.ReadInputFiles(args.treeFile, args.geneTreeFile, args.gtFormat)
	if err != nil {
		log.Fatalf("%s %s\n", ErrMessage, err)
	}
	switch args.command {
	case Infer:
		log.Println("running infer...")
		td, results, err := infer.CAMUS(tre, geneTrees.Trees)
		if err != nil {
			log.Fatalf("%s %s\n", ErrMessage, err)
		}
		for _, branches := range results {
			fmt.Println(gr.MakeNetwork(td, branches).Newick())
		}
	case Score:
		log.Println("running score...")
		network, err := pr.ConvertToNetwork(tre)
		if err != nil {
			log.Fatalf("%s %s\n", ErrMessage, err)
		}
		scores, err := score.CalculateReticulationScore(network, geneTrees.Trees)
		if err != nil {
			log.Fatalf("%s %s\n", ErrMessage, err)
		}
		pr.WriteBranchScoresToCSV(scores, geneTrees.Names)
	default:
		panic(fmt.Sprintf("invalid command (%d)", args.command))
	}
}
