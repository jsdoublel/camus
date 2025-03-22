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

There is three methods for installing go. 

### Using Go's package manager

After installing Go, run

```sh
go install github.com/jsdoublel/camus@latest
```

this will install the CAMUS executable to `$GOPATH/bin`. You can check where
that is using the command `go env GOPATH`; if you want to run CAMUS from any
directory make sure you add `$GOPATH/bin` to your `PATH` environmental variable.

### Build from source

Once Go is installed, clone the repository, download dependencies, and build
the project with the following commands

```sh
git clone https://github.com/jsdoublel/camus.git
cd camus
go mod download
go build
```

### Download the binaries

The GitHub [releases](https://github.com/jsdoublel/camus/releases) contain
binaries compiled on my personal computer (i.e., Linux x86_64). 

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

