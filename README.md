# CAMUS

[![Go Reference](https://pkg.go.dev/badge/github.com/jsdoublel/camus.svg)](https://pkg.go.dev/github.com/jsdoublel/camus)
[![GitHub go.mod Go version](https://img.shields.io/github/go-mod/go-version/jsdoublel/camus?logo=go)](https://go.dev/)
[![Go Report Card](https://goreportcard.com/badge/github.com/jsdoublel/camus)](https://goreportcard.com/report/github.com/jsdoublel/camus)
[![build](https://github.com/jsdoublel/camus/actions/workflows/go.yml/badge.svg)](https://github.com/jsdoublel/camus/actions)
[![GitHub License](https://img.shields.io/github/license/jsdoublel/camus)](https://github.com/jsdoublel/camus/blob/main/LICENSE.txt)

CAMUS (Constrained Algorithm Maximizing qUartetS) is a dynamic programming
algorithm for inferring level-1 phylogenetic networks from quartets and a
constraint tree.

## Algorithm

Given a rooted, binary, constraint tree $T$ and a list of input trees
$\mathcal{G}$, the CAMUS algorithm finds the level-1 networks inducing the
maximum number of quartets from $\mathcal{G}$ that contains the constraint tree
$T$

CAMUS finds networks for different values of $k$, where $k$ is the number of
edges added to the network. For example, if the best possible network that
contains $T$ has $m$ edges, then CAMUS will output $m$ networks $N_1, N_2,
\cdots, N_m$, where $N_i$ is optimal under the constraint that it contains
exactly $i$ edges.

CAMUS has the following inputs and outputs:

**Input**

- *Constraint Tree:* Rooted, binary, tree in newick format without duplicate
  labels.
- *Gene Trees:* List of trees in newick format, containing only labels from the
  constraint tree.

**Output**

- *Output Network:* Level-1 networks written in extended newick format.

CAMUS  should be invoked with the constraint tree file path and gene trees file 
path as positional arguments in that order; the output network and logging 
information is written to files with a prefix that can optionally be set 
with the `-o` flag.

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

Once Go is installed, clone the repository and build the project with the 
following commands:

```
git clone https://github.com/jsdoublel/camus.git
cd camus
go build
```

### Download the binaries

The GitHub [releases](https://github.com/jsdoublel/camus/releases) contain
binaries compiled for Windows, macOS, and Linux (for both arm64 and amd64).

## Usage

```
camus [ -f <format> | -o <output> | -t <threshold> | -n <threads> | -h | -v | ... ] <const_tree> <gene_trees>
```

There are two positional arguments indicating the inputs. Additionally, there
are the following flags that precedes the positional arguments.

**Basic Flags**

- `-f format [ newick | nexus ] (default "newick")` sets the format of the input gene tree file
- `-t threshold [0, 1] (default 0.5)` quartet filtering threshold
- `-n num_procs` number of parallel processes
- `-o prefix` output prefix
- `-h` prints usage information and exits
- `-hh` prints extended usage information and exits
- `-v` prints software version and exits

**Experimental Flags**
 
- `-sm mode [ max | norm | sym ] (default "max")` sets the score mode
- `-a alpha` parameter that adjusts penalty in ``sym" score mode
- `-asSet` quartet count is calculated as a set (counts total unique quartet topologies)
- `-s threshold` collapse edges in gene trees with support less than threshold value
- `-q mode [0, 2] (default 0)` quartet filtering mode
  
### Quartet Filter Mode

Quartet filtering mode filters out less frequent quartet topologies. Mode `-q
0` disables quartet filtering; `-q 1` applies a less restrictive quartet
filtering, and `-q 2` is the most restrictive and recommended quartet
filtering.

### Score Modes

Score Modes are various modifications to the optimization score beyond 
simple maximization. These are experimental and maximization is recommended.



