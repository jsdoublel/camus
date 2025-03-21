# CAMUS 

CAMUS (Constrained Algorithm Maximizing qUartetS) is a dynamic programming
algorithm for inferring level-1 phylogenetic networks from quartets and a
constraint tree.

Specifically, given a rooted, binary, constraint tree $T$ and a list of input
trees $\mathcal{G}$, CAMUS finds the level-1 network inducing the maximum
number of quartets from $\mathcal{G}$ that respects the constraint tree $T$.

## Installation

CAMUS should be able to build on any operating system, though it has only been
tested on Linux. It is implemented in the Go programming language (which can be
installed [here](https://go.dev/doc/install)) and additionally requires
[gotree](https://github.com/evolbioinfo/gotree).

Once Go is installed, further dependencies can be installed with the command
`go mod download` from the project directory. Then simply run `go build` to
build the project.

## Usage

**Input**

- *Constraint Tree:* Rooted, binary, tree in newick format without duplicate
  labels.
- *Gene Trees:* List of trees in newick format, containing only labels from the
  constraint tree.

**Output**

- *Output Network:* Level-1 network written in extended newick format.

CAMUS is invoked with the constraint tree file path and gene trees file path as
positional arguments in that order; the output network  and logging information
is written to `stdout` and `stderr` respectivly.

Here is an example of how one might run CAMUS using example data from this repository.

```sh
camus testdata/benchmark/constraint.nwk testdata/benchmark/gene-trees.nwk > out.nwk 2> log.txt
```

