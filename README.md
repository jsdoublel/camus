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

## TODO

- [x] code algorithm
- [ ] fix error message when gene tree has multiple of the same label
- [x] test case when edge connects down to ancestor 
- [ ] test case when two edges connect to the same edge in the constraint tree
- [ ] fix nil pointer if `isBinary` probably because node had only one child
- [x] improve logging
