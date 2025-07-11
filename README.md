# CAMUS

[![build](https://github.com/jsdoublel/camus/actions/workflows/go.yml/badge.svg)](https://github.com/jsdoublel/camus/actions)
[![Go Reference](https://pkg.go.dev/badge/github.com/jsdoublel/camus.svg)](https://pkg.go.dev/github.com/jsdoublel/camus)
[![GitHub go.mod Go version](https://img.shields.io/github/go-mod/go-version/jsdoublel/camus?logo=go)](https://go.dev/)
[![GitHub License](https://img.shields.io/github/license/jsdoublel/camus)](https://github.com/jsdoublel/camus/blob/main/LICENSE.txt)

CAMUS (Constrained Algorithm Maximizing qUartetS) is a dynamic programming
algorithm for inferring level-1 phylogenetic networks from quartets and a
constraint tree.

## Installation

CAMUS should be able to build on any operating system, though it has only been
tested on Linux. It is implemented in the Go programming language (which can be
installed [here](https://go.dev/doc/install)).

There are three methods for installing CAMUS. 

### Using Go's package manager

After installing Go, run:

```
go install github.com/jsdoublel/camus@latest
```

This will install the CAMUS executable to `$GOPATH/bin`. You can check where
that is using the command `go env GOPATH`; if you want to run CAMUS from any
directory make sure you add `$GOPATH/bin` to your `PATH` environmental variable.

### Build from source

Once Go is installed, clone the repository and build the project with the following commands:

```
git clone https://github.com/jsdoublel/camus.git
cd camus
go build
```

### Download the binaries

The GitHub [releases](https://github.com/jsdoublel/camus/releases) contain
binaries compiled on my personal computer (i.e., Linux x64). 

## Usage

```
camus [ -f <format> | -h | -v ] <command> <tree> <gene_trees>
```

There are two main commands for CAMUS: `infer` and `score`, followed by two positional arguments 
indicating the inputs. Addtionally, there are the following flags that preceed the positional arguments.

- `-f format [ newick | nexus ] (default "newick")` sets the format of the input gene tree file
- `-q quartet filter mode [0, 2] (default 0)` defines quartet filtering mode
- `-h` prints usage information and exits
- `-v` prints software version and exits
  
#### Quartet filter mode

Quartet filter mode is a mode for filtering out lower frequency quartet topologies. Since there are three quartet topologies for every for taxa, valid values range from 0 to 2. Mode 0 filters out no quartets, mode 1 filters out the single least fequently seen quartet, and mode 2 filters out the two least fequently seen quartets. 

> Quartet filter mode currently only works with the `infer` command.

### Infer

```
camus infer <constraint_tree> <gene_trees>
```

Given a rooted, binary, constraint tree $T$ and a list of input trees
$\mathcal{G}$, the CAMUS algorithm finds the level-1 networks inducing the
maximum number of quartets from $\mathcal{G}$ that contains the constraint tree
$T$

CAMUS finds networks for different values of $k$, where $k$ is the number of
edges added to the network. For example, if the best possible network that
contains $T$ has $m$ edges, then CAMUS will output $m$ networks $N_1, N_2,
\cdots, N_m$, where $N_i$ is optimal under the constraint that it contains
exactly $i$ edges.

Since these added edges are directed edges, it is a requirement that any
quartet that contributes to the maximum contain exactly one taxon from the clade
below by $w$—where $w$ is the vertex that the new edge points towards.

CAMUS has the following inputs and outputs:

**Input**

- *Constraint Tree:* Rooted, binary, tree in newick format without duplicate
  labels.
- *Gene Trees:* List of trees in newick format, containing only labels from the
  constraint tree.

**Output**

- *Output Network:* Level-1 networks written in extended newick format.

CAMUS, when run with the `infer` command, should be invoked with the constraint
tree file path and gene trees file path as positional arguments in that order;
the output network  and logging information is written to `stdout` and `stderr`
respectively.

Here is an example of how one might run CAMUS using example data from this
repository:

```
camus infer testdata/large/constraint.nwk testdata/large/gene-trees.nwk > out.nwk 2> log.txt
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
we mean the set of quartets where each taxon $t$ in the quartet corresponds to a
unique node in $\gamma_r$ that can be reached by a path from $t$ to that unique
node, without passing through any other vertices in $\gamma_r$. 

We also define $\mathcal{Q}_{\gamma_r}'$ where $\forall q' \in \mathcal{Q}\_{\gamma_r}', 
\exists q \in \mathcal{Q}\_{\gamma_r}$ where $\mathcal{L}(q') = \mathcal{L}(q)$, where 
$\mathcal{L}(q)$ is the set of taxa in $q$ (in short, $\mathcal{Q}\_{\gamma_r}'$ is 
the set of all possible quartets on the same sets of taxa as the quartets in 
$\mathcal{Q}\_{\gamma_r}$).

Thus, given a quartet set $\mathcal{Q}\_g$ from a gene tree $g$, we can
calculate the support of $g$ for $r$ as 

$$S_{g,r} = \frac{|(\mathcal{Q}\_{\gamma_r} \cap
\mathcal{Q}\_g) \setminus \mathcal{Q}\_T|}{|\mathcal{Q}\_{\gamma_r}' \cap \mathcal{Q}\_g|}$$

where $\mathcal{Q}\_T$ is the set of quartets in the constraint tree.
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

As with `infer`, `score` uses positional arguments for the inputs and outputs 
the result and logging information to `stdout` and `stderr` respectively. 

Here is an example of how one might run CAMUS using example data from this
repository:

```
camus score testdata/large/network.nwk testdata/large/gene-trees.nwk > out.csv 2> log.txt
```
