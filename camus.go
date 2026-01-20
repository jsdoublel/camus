/*
CAMUS (Constrained Algorithm Maximizing qUartetS) is a dynamic programming
algorithm for inferring level-1 phylogenetic networks from quartets and a
constraint tree.

# MIT License

# Copyright (c) 2026 James Willson

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

usage: camus [flags]... <const_tree_file> <gene_tree_file>

positional arguments:

	<tree_file>			constraint newick tree
	<gene_tree_file>	gene tree newick file

flags:

	-f format
	  	gene tree format [newick|nexus] (default "newick")
	-h	prints short help and exits
	-hh
	  	prints help with experimental features and exits
	-n int
	  	number of parallel processes
	-o string
	  	output prefix
	-t float
	  	threshold for quartet filter [0, 1] (default 0.5)
	-v	prints version number and exits

examples:

	camus -o output-name constraint.nwk gene-trees.nwk
*/
package main

import (
	"bytes"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"slices"
	"strings"
	"time"

	gr "github.com/jsdoublel/camus/internal/graphs"
	in "github.com/jsdoublel/camus/internal/infer"
	pr "github.com/jsdoublel/camus/internal/prep"
	sc "github.com/jsdoublel/camus/internal/score"
)

const (
	Version      = "v1.0.0"
	ErrorMessage = "camus encountered an error ::"
	TimeFormat   = "2006-01-02_15-04-05"

	DefaultFormat     = "newick"
	DefaultScoreMode  = "max"
	DefaultQMode      = 2
	DefaultMinSupport = 0
	DefaultThreshold  = 0.5
	DefaultAlpha      = 0.1
)

var experimentalFlags = []string{"a", "asSet", "q", "s", "sm"}

type Args struct {
	prefix       string          // output prefix
	gtFormat     pr.Format       // gene tree file format
	treeFile     string          // constraint or network tree file
	geneTreeFile string          // gene trees
	inferOpts    in.InferOptions // camus options
}

func Usage(extended bool) {
	fmt.Fprint(flag.CommandLine.Output(), // nolint
		"usage: camus [flags]... <const_tree_file> <gene_tree_file>\n",
		"\n",
		"positional arguments:\n\n",
		"  <tree_file>\t\tconstraint newick tree\n",
		"  <gene_tree_file>\tgene tree newick file\n",
		"\n",
		"flags:\n\n",
	)
	if extended {
		flag.PrintDefaults()
	} else {
		shortFlags := flag.NewFlagSet("short", flag.ContinueOnError)
		flag.VisitAll(func(f *flag.Flag) {
			if !slices.Contains(experimentalFlags, f.Name) {
				shortFlags.Var(f.Value, f.Name, f.Usage)
				shortFlags.Lookup(f.Name).DefValue = f.DefValue
			}
		})
		shortFlags.PrintDefaults()
	}
	fmt.Fprint(flag.CommandLine.Output(), // nolint
		"\n",
		"examples:\n\n",
		"\tcamus -o output-name constraint.nwk gene-trees.nwk\n\n",
	)
}

func parseArgs() Args {
	flag.Usage = func() {
		Usage(false)
	}
	format, ok := pr.ParseFormat[DefaultFormat]
	if !ok {
		panic(fmt.Sprintf("bad default format %s", DefaultFormat))
	}
	flag.Var(&format, "f", "gene tree `format` [newick|nexus] (default \"newick\")")
	prefix := flag.String("o", "", "output prefix")
	scoreMode := flag.String("sm", DefaultScoreMode, "score `mode` [max|norm|sym]")
	mode := flag.Int("q", DefaultQMode, "quartet filter mode number [0, 2]")
	supp := flag.Float64("s", DefaultMinSupport, "collapse edges in gene trees with support less than value (default 0)")
	thresh := flag.Float64("t", DefaultThreshold, "threshold for quartet filter [0, 1]")
	alpha := flag.Float64("a", DefaultAlpha, "parameter to adjust penalty for \"sym\" score mode, from (0, 1]")
	asSet := flag.Bool("asSet", false, "quartet count is calculated as a set (one point per unique topology)")
	help := flag.Bool("h", false, "prints short help and exits")
	hhelp := flag.Bool("hh", false, "prints help with experimental features and exits")
	ver := flag.Bool("v", false, "prints version number and exits")
	nprocs := flag.Int("n", 0, "number of parallel processes")
	flag.Parse()
	if *help {
		Usage(false)
		os.Exit(0)
	}
	if *hhelp {
		Usage(true)
		os.Exit(0)
	}
	if *ver {
		fmt.Printf("camus %s\n", Version)
		os.Exit(0)
	}
	if flag.NArg() != 2 {
		parserError("two positional arguments required: <const_tree> <gene_tree_file>")
	}
	scorer, ok := sc.ParseScorer[*scoreMode]
	if !ok {
		parserError(fmt.Sprintf("\"%s\" is not a valid score mode: valid score modes are \"max\", \"norm\", and \"sym\"", *scoreMode))
	}
	qOpts, err := pr.SetQuartetFilterOptions(*mode, *thresh)
	if err != nil {
		parserError(err.Error())
	}
	inferOpts, err := in.MakeInferOptions(*nprocs, qOpts, *supp, scorer, *asSet, *alpha)
	if err != nil {
		parserError(err.Error())
	}
	return Args{
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
	Usage(false)
	os.Exit(1)
}

func defaultPrefix() string {
	parseName := func(s string) string {
		parts := strings.Split(s, string(os.PathSeparator))
		parts = strings.Split(parts[len(parts)-1], ".")
		if len(parts) > 1 {
			return strings.Join(parts[:len(parts)-1], ".")
		}
		return parts[0]
	}
	inputs := fmt.Sprintf("%s_%s", parseName(flag.Arg(0)), parseName(flag.Arg(1)))
	return fmt.Sprintf("camus_%s_%s", inputs, time.Now().Local().Format(TimeFormat))
}

func main() {
	var exit int
	defer func() {
		os.Exit(exit)
	}()
	buf := &bytes.Buffer{} // capture pre logfile setup logging
	log.SetFlags(log.LstdFlags | log.Lmicroseconds)
	log.SetOutput(io.MultiWriter(os.Stderr, buf))
	args := parseArgs()
	if args.prefix == "" {
		args.prefix = defaultPrefix()
		log.Printf("output prefix was not set, using \"%s\"", args.prefix)
	}
	if logf, err := os.Create(fmt.Sprintf("%s.log", args.prefix)); err == nil {
		logf.Write(buf.Bytes()) // nolint
		log.SetOutput(io.MultiWriter(os.Stderr, logf))
		defer func() {
			log.SetOutput(os.Stderr)
			_ = logf.Close()
		}()
	} else {
		log.Printf("failed to create log file %s.log, %s", args.prefix, err) // should continue to log to stderr
	}
	log.Printf("camus %s", Version)
	log.Printf("invoked as: camus %s", strings.Join(os.Args[1:], " "))
	if err := run(args); err != nil {
		log.Printf("%s %s", ErrorMessage, err)
		exit = 1
	}
}

func run(args Args) error {
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
