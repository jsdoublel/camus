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
  -g string
        gene tree file
  -o string
        output extended newick file
  -t string
        constraint tree
```

## TODO

- [x] code algorithm
- [x] generally improve error handling
- [x] test case when edge connects down to ancestor 
- [x] test case when two edges connect to the same edge in the constraint tree
- [ ] fix nil pointer if `isBinary` probably because node had only one child (try to replicate)
- [ ] fix error message when gene tree has multiple of the same label or if
  there is a label in the gene tree not in the constraint tree
- [x] improve logging
- [x] add benchmark
- [ ] update comments for functions to at least make sure they're accurate
- [ ] check empty inputs
- [ ] finish up readme
- [ ] more than one constraint tree
- [ ] improper newick format (both gene tree and constraint tree)
- [ ] leaves appear in the cosntraint tree, but not in any input trees
- [ ] leaves appear in input trees, but not in constraint tree
- [ ] profile cpu and memory usage
- [x] look into creating proper module url
