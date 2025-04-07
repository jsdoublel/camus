/*
CAMUS (Constrained Algorithm Maximizing qUartetS) is a dynamic programming
algorithm for inferring level-1 phylogenetic networks from quartets and a
constraint tree.

usage: camus [-h| -v | -f <format> ] <command> <tree> <gene_trees>

commands:

	infer		find level-1 network given constraint tree and gene trees
	score		score each reticulation branch with respects to gene trees

positional arguments:

	<tree>			constraint newick tree (infer) or network (score)
	<gene_trees>	gene tree newick file

flags:

	-f string
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

	"github.com/mattn/go-isatty"

	"github.com/jsdoublel/camus/graphs"
	"github.com/jsdoublel/camus/infer"
	"github.com/jsdoublel/camus/prep"
	"github.com/jsdoublel/camus/score"
)

var version = "v0.2.3"

type args struct {
	command      string // infer or score
	treeFile     string // constraint or network tree file
	geneTreeFile string // gene trees
	gtFormat     string // gene tree file format
}

var ErrMessage = red("error:")

func parseArgs() args {
	flag.Usage = func() {
		fmt.Fprint(os.Stderr,
			"usage: camus [-h| -v | -f <format> ] <command> <tree> <gene_trees>\n",
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
	help := flag.Bool("h", false, "prints this message and exits")
	ver := flag.Bool("v", false, "prints version number and exits")
	format := flag.String("f", "newick", "gene tree format [ newick | nexus ]")
	flag.Parse()
	if *help {
		flag.Usage()
		os.Exit(0)
	}
	if *ver {
		fmt.Printf("CAMUS version %s", version)
		os.Exit(0)
	}
	if *format != "nexus" && *format != "newick" {
		parserError(fmt.Sprintf("%s \"%s\" is not a valid gene tree file format", ErrMessage, *format))
	}
	if flag.Arg(0) == "" {
		parserError(fmt.Sprintf("%s no command given: either infer or score required", ErrMessage))
	}
	if flag.Arg(0) != "infer" && flag.Arg(0) != "score" {
		parserError(fmt.Sprintf("%s \"%s\" is not a valid command", ErrMessage, flag.Arg(0)))
	}
	if flag.NArg() != 3 {
		parserError(fmt.Sprintf("%s two positional arguments required: <tree> <gene_tree_file>", ErrMessage))
	}
	return args{
		command:      flag.Arg(0),
		treeFile:     flag.Arg(1),
		geneTreeFile: flag.Arg(2),
		gtFormat:     *format,
	}
}

// prints message, usage, and exits (statis code 1)
func parserError(message string) {
	fmt.Fprintln(os.Stderr, message)
	flag.Usage()
	os.Exit(1)
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
	tre, geneTrees, err := prep.ReadInputFiles(args.treeFile, args.geneTreeFile, args.gtFormat)
	if err != nil {
		log.Fatalf("%s %s\n", ErrMessage, err)
	}
	switch args.command {
	case "infer":
		td, branches, err := infer.CAMUS(tre, geneTrees.Trees)
		if err != nil {
			log.Fatalf("%s %s\n", ErrMessage, err)
		}
		fmt.Println(graphs.MakeNetwork(td, branches).Newick())
	case "score":
		network, err := prep.ConvertToNetwork(tre)
		if err != nil {
			log.Fatalf("%s %s\n", ErrMessage, err)
		}
		scores, err := score.CalculateReticulationScore(network, geneTrees.Trees)
		if err != nil {
			log.Fatalf("%s %s\n", ErrMessage, err)
		}
		prep.WriteBranchScoresToCSV(scores, geneTrees.Names)
	default:
		panic(fmt.Sprintf("allowed invalid command %s", args.command))
	}
}
