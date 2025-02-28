# CAMUS 

CAMUS (Constrained Algorithm Maximizing qUartetS) is a dynamic programming algorithm for inferring level-1 phylogenetic networks from quartets and a constraint tree.

## Installation

Requires [gotree](https://github.com/evolbioinfo/gotree). Dependencies can be installed with the command `go mod download` from the project directory. Then simply run `go build` to build the project.

## Usage

```
  -g string
        gene tree file
  -o string
        output extended newick file
  -t string
        constraint tree
```
