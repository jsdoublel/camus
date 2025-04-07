# CAMUS
[![build](https://github.com/jsdoublel/camus/actions/workflows/go.yml/badge.svg)](https://github.com/jsdoublel/camus/actions)

CAMUS (Constrained Algorithm Maximizing qUartetS) is a dynamic programming
algorithm for inferring level-1 phylogenetic networks from quartets and a
constraint tree.


## Installation

CAMUS should be able to build on any operating system, though it has only been
tested on Linux. It is implemented in the Go programming language (which can be
installed [here](https://go.dev/doc/install)).

Additionally, CAMUS depends on the following three go modules:

- [bitset](https://github.com/bits-and-blooms/bitset)
- [go-isatty](https://github.com/mattn/go-isatty).
- [gotree](https://github.com/evolbioinfo/gotree) 

There are three methods for installing CAMUS. 

### Using Go's package manager

After installing Go, run

```
go install github.com/jsdoublel/camus@latest
```

this will install the CAMUS executable to `$GOPATH/bin`. You can check where
that is using the command `go env GOPATH`; if you want to run CAMUS from any
directory make sure you add `$GOPATH/bin` to your `PATH` environmental variable.

### Build from source

Once Go is installed, clone the repository, download dependencies, and build
the project with the following commands

```
git clone https://github.com/jsdoublel/camus.git
cd camus
go mod download
go build
```

### Download the binaries

The GitHub [releases](https://github.com/jsdoublel/camus/releases) contain
binaries compiled on my personal computer (i.e., Linux x86_64). 

## Usage

```
camus [-h| -v | -f <format> ] <command> <tree> <gene_trees>
```

There are two main commands for CAMUS: `infer` and `score`, followed by two positional arguments indicating the inputs. Addtionally, there are the following flags that preceed the positional arguments.

- `-f string [ newick | nexus ] (default "newick")` sets the format of the input gene tree file
- `-h` prints usage information and exits
- `-v` prints software version and exits

### Infer

```
camus infer <constraint_tree> <gene_trees>
```

Given a rooted, binary, constraint tree $T$ and a list of input trees
$\mathcal{G}$, the CAMUS algorithm finds the level-1 network inducing the
maximum number of quartets from $\mathcal{G}$ that respects the constraint tree
$T$.

Since these added edges are directed edges, it is a requirement that any
quartet that contributes to the maximum contain exactly one taxa from the clade
below by $w$—where $w$ is the vertex that the new edge points towards.

CAMUS has the following inputs and outputs:

**Input**

- *Constraint Tree:* Rooted, binary, tree in newick format without duplicate
  labels.
- *Gene Trees:* List of trees in newick format, containing only labels from the
  constraint tree.

**Output**

- *Output Network:* Level-1 network written in extended newick format.

CAMUS, when run with the `infer` command, should be invoked with the constraint
tree file path and gene trees file path as positional arguments in that order;
the output network  and logging information is written to `stdout` and `stderr`
respectively.

Here is an example of how one might run CAMUS using example data from this
repository:

```
camus infer testdata/benchmark/constraint.nwk testdata/benchmark/gene-trees.nwk > out.nwk 2> log.txt
```

### Score

```
camus score <network> <gene_trees>
```

When using the `score` command, CAMUS returns a CSV file indicating the support
for each branch for each gene. Support is calculated as the ratio of the
relevant quartets from the gene that support the branch.

Specifically, for each new reticulation edge $r$ creating a cycle $\gamma_r$,
we have $\mathcal{Q}\_{\gamma_r}$—the set of quartets "in" the cycle. By this
we mean the set of quartets where each taxa $t$ in the quartet corresponds to a
unique node in $\gamma_r$ that can be reached by a path from $t$ to that unique
node, without passing through any other vertices in $\gamma_r$. 

We also define $\mathcal{Q}_{\gamma_r}'$ where $\forall q' \in \mathcal{Q}\_{\gamma_r}', \exists q \in \mathcal{Q}\_{\gamma_r}$ where $\mathcal{L}(q') = \mathcal{L}(q)$, where $\mathcal{L}(q)$ is the set of taxa in $q$ (in short, $\mathcal{Q}\_{\gamma_r}'$ is the set of quartets that are relevant for $\gamma_t$).

Thus, given a quartet set $\mathcal{Q}\_g$ from a gene tree $g$, we can
calculate the support of $g$ for $r$ as 

$$S_{g,r} = \frac{|\mathcal{Q}\_{\gamma_r} \cap
\mathcal{Q}\_g|}{|\mathcal{Q}\_{\gamma_r}' \cap \mathcal{Q}\_g|}$$

As can be seen, if the inputted gene trees were simply quartets (four leaf
newick trees), the only valid results would be 0, 1, or NaN.

The `score` command for CAMUS has the following inputs and outputs:

**Input**

- *Network:* Rooted, fully-resolved, network in extended newick format.
- *Gene Trees:* List of trees in newick format, containing only labels from the
  constraint tree.

**Output**

- *Scores:* CSV file where the columns are the reticulation branch labels and
  the rows correspond to different genes (labeled by the line number they
  appear on in the input).

As with `infer`, `score` uses positional arguments for the inputs and outputs the result and logging information to `stdout` and `stderr` respectively. 

Here is an example of how one might run CAMUS using example data from this
repository:

```
camus score testdata/benchmark/network.nwk testdata/benchmark/gene-trees.nwk > out.csv 2> log.txt
```
