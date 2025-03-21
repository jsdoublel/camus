# CAMUS 

CAMUS (Constrained Algorithm Maximizing qUartetS) is a dynamic programming
algorithm for inferring level-1 phylogenetic networks from quartets and a
constraint tree.

## Installation

Requires [gotree](https://github.com/evolbioinfo/gotree). Dependencies can be
installed with the command `go mod download` from the project directory. Then
simply run `go build` to build the project.

## Usage

```
usage: camus [-h] <constraint_tree> <gene_tree>

positional arguments (required):
  <constraint_tree>        constraint newick tree
  <gene_tree>              gene tree newick file

flags:
  -h	prints this message and exits

example:
  github.com/jsdoublel/camus contraint.nwk gene-trees.nwk > out.nwk 2> log.txt
```

