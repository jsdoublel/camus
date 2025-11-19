package prep

import (
	"bufio"
	"bytes"
	"encoding/csv"
	"errors"
	"fmt"
	"image/color"
	"io"
	"log"
	"math"
	"os"
	"slices"
	"strconv"
	"strings"

	gr "github.com/jsdoublel/camus/internal/graphs"

	"github.com/evolbioinfo/gotree/io/newick"
	"github.com/evolbioinfo/gotree/io/nexus"
	"github.com/evolbioinfo/gotree/tree"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/vg"
	"gonum.org/v1/plot/vg/draw"
)

var (
	ErrInvalidFile     = errors.New("invalid file")
	ErrInvalidFormat   = errors.New("invalid format")
	ErrNoReticulations = errors.New("no reticulations")
	ErrWritingFile     = errors.New("error writing file")

	plotLineColor  = color.RGBA{R: 37, G: 150, B: 190, A: 255}
	plotMarkerShap = draw.SquareGlyph{}
)

type Format int

const (
	Newick Format = iota
	Nexus

	plotH = 4 * vg.Inch
	plotW = 6 * vg.Inch

	maxTicks = 10
)

var ParseFormat = map[string]Format{
	"newick": Newick,
	"nexus":  Nexus,
}

func (f *Format) Set(s string) error {
	if format, ok := ParseFormat[s]; ok {
		*f = format
		return nil
	}
	return fmt.Errorf("\"%s\" is not a valid gene tree file format", s)
}

func (f Format) String() string {
	for s, fr := range ParseFormat {
		if fr == f {
			return s
		}
	}
	panic(fmt.Sprintf("format (%d) does not exist", f))
}

type GeneTrees struct {
	Trees []*tree.Tree // gene trees
	Names []string     // gene names
}

// Reads in and validates constraint tree and gene tree input files.
// Returns an error if the newick format is invalid, or the file is invalid for
// some other reason (e.g., more than one constraint tree)
func ReadInputFiles(treeFile, genetreesFile string, format Format) (*tree.Tree, *GeneTrees, error) {
	flags := log.Flags()
	lout := log.Writer()
	log.SetOutput(io.Discard) // don't log this bit as gotree can be noisy and lead to thousands of log messages
	defer func() {
		log.SetOutput(lout)
		log.SetFlags(flags)
	}()
	tre, err := readTreeFile(treeFile)
	if err != nil {
		return nil, nil, err
	}
	genetrees, err := readGeneTreesFile(genetreesFile, format)
	if err != nil {
		return nil, nil, err
	}
	return tre, genetrees, nil
}

// reads and validates constraint tree file
func readTreeFile(treeFile string) (*tree.Tree, error) {
	treBytes, err := os.ReadFile(treeFile)
	if err != nil {
		return nil, fmt.Errorf("error reading tree file: %w", err)
	}
	treBytes = bytes.TrimSpace(treBytes)
	if bytes.Count(treBytes, []byte{byte('\n')}) != 0 || len(treBytes) == 0 {
		return nil, fmt.Errorf("%w, there should only be exactly one newick tree in tree file %s",
			ErrInvalidFile, treeFile)
	}
	tre, err := newick.NewParser(bytes.NewReader(treBytes)).Parse()
	if err != nil {
		return nil, fmt.Errorf("%w, error parsing tree newick string from %s: %s",
			ErrInvalidFormat, treeFile, err.Error())
	}
	tre.ClearLengths(true, true)
	tre.ClearComments()
	tre.ClearSupports()
	return tre, nil
}

// reads and validates gene tree file
func readGeneTreesFile(genetreesFile string, format Format) (*GeneTrees, error) {
	file, err := os.Open(genetreesFile)
	if err != nil {
		return nil, fmt.Errorf("error opening %s, %w", genetreesFile, err)
	}
	defer func() {
		if err := file.Close(); err != nil {
			panic(fmt.Sprintf("could not close file %s, %s", genetreesFile, err))
		}
	}()
	geneTreeList := make([]*tree.Tree, 0)
	geneTreeNames := make([]string, 0)
	switch format {
	case Newick:
		scanner := bufio.NewScanner(file)
		for i := 0; scanner.Scan(); i++ {
			line := bytes.TrimSpace(scanner.Bytes())
			if line != nil {
				genetree, err := newick.NewParser(bytes.NewReader(line)).Parse()
				if err != nil {
					return nil, fmt.Errorf("%w, error reading gene tree on line %d in %s: %s",
						ErrInvalidFormat, i, genetreesFile, err.Error())
				}
				geneTreeList = append(geneTreeList, genetree)
			}
		}
		if len(geneTreeList) < 1 {
			return nil, fmt.Errorf("%w, empty gene tree file %s", ErrInvalidFile, genetreesFile)
		}
		geneTreeNames = make([]string, 0)
		for i := range len(geneTreeList) {
			geneTreeNames = append(geneTreeNames, strconv.Itoa(i+1))
		}
	case Nexus:
		nex, err := nexus.NewParser(file).Parse()
		if err != nil {
			return nil, fmt.Errorf("%w, error reading gene tree nexus file %s: %s",
				ErrInvalidFormat, genetreesFile, err.Error())
		}
		nex.IterateTrees(func(s string, t *tree.Tree) {
			geneTreeList = append(geneTreeList, t)
			geneTreeNames = append(geneTreeNames, s)
		})
	default:
		return nil, fmt.Errorf("%w, not a valid file format", ErrInvalidFile)
	}
	return &GeneTrees{Trees: geneTreeList, Names: geneTreeNames}, nil
}

// Read in extended newick file and make network
func ConvertToNetwork(ntw *tree.Tree) (network *gr.Network, err error) {
	if !ntw.Rooted() {
		return nil, fmt.Errorf("network is %w", ErrUnrooted)
	}
	if !NetworkIsBinary(ntw) {
		return nil, fmt.Errorf("network is %w", ErrNonBinary)
	}
	ret := make(map[string]gr.Branch)
	var errNode *tree.Node
	ntw.PostOrder(func(cur, prev *tree.Node, e *tree.Edge) (keep bool) {
		if errNode != nil {
			return true
		}
		if strings.Contains(cur.Name(), "#") {
			branch := ret[cur.Name()]
			var v *tree.Node
			if cur.Tip() {
				p, err := prev.Parent()
				if err != nil && err.Error() != "The node has no parent : May be the root?" {
					panic(fmt.Sprintf("err from backbone tree: %s", err))
				}
				for _, n := range prev.Neigh() {
					if n != cur && n != p {
						v = n
					}
				}
				if branch.IDs[gr.Ui] != 0 || v == nil {
					errNode = cur
					return true
				}
				branch.IDs[gr.Ui] = v.Id()
			} else {
				for _, n := range cur.Neigh() {
					if n != cur && n != prev {
						v = n
					}
				}
				if branch.IDs[gr.Wi] != 0 || v == nil {
					errNode = cur
					return true
				}
				branch.IDs[gr.Wi] = v.Id()
			}
			ret[cur.Name()] = branch
		}
		return true
	})
	if errNode != nil {
		return nil, fmt.Errorf("%w, too many or invalid matching reticulation label %s",
			ErrInvalidFormat, errNode.Name())
	}
	if len(ret) == 0 {
		return nil, fmt.Errorf("%w - not a network", ErrNoReticulations)
	}
	for label, branch := range ret {
		if branch.IDs[gr.Ui] == 0 || branch.IDs[gr.Wi] == 0 { // assumes root node is not labeled as reticulation
			return nil, fmt.Errorf("%w, label %s is unmatched", ErrInvalidFormat, label)
		}
	}
	if err := ntw.UpdateTipIndex(); err != nil {
		return nil, fmt.Errorf("network %w", ErrMulTree)
	}
	return &gr.Network{NetTree: ntw, Reticulations: ret}, nil
}

// Write DP results csv file to writer.
//
// There are three columns: "Number of Branches", "Quartet Satisfied Percent", "Extended Newick"
func WriteDPResultsToCSV(td *gr.TreeData, newicks []string, qsat []float64, w io.Writer) (err error) {
	if len(newicks) != len(qsat) {
		panic(fmt.Sprintf("there should be a set of branches for every optimal score, %+v %+v", newicks, qsat))
	}
	data := make([][]string, len(newicks)+2)
	data[0] = []string{"Number of Branches", "Quartet Satisfied Percent", "Extended Newick"}
	data[1] = []string{strconv.FormatInt(0, 10), strconv.FormatFloat(0, 'f', -1, 64), td.Newick()}
	for i := range len(newicks) {
		data[i+2] = []string{
			strconv.FormatInt(int64(i+1), 10),
			strconv.FormatFloat(qsat[i], 'f', -1, 64),
			newicks[i],
		}
	}
	writer := csv.NewWriter(w)
	defer func() {
		writer.Flush()
		if err == nil {
			err = writer.Error()
		} else if writer.Error() != nil {
			log.Printf("error when flushing output csv, %s", writer.Error())
		}
	}()
	if err = writer.WriteAll(data); err != nil {
		err = fmt.Errorf("%w, %s", ErrWritingFile, err)
		return
	}
	return
}

func WriteResultsLineplot(qstat []float64, prefix string) error {
	p := plot.New()
	p.X.Label.Text = "Number of Reticulations"
	p.Y.Label.Text = "Percent of Quartets Not Satisfied"
	p.X.Min = 0
	p.X.Max = float64(len(qstat))
	p.X.Tick.Marker = plot.TickerFunc(func(_, max float64) []plot.Tick {
		step := 1
		if int(max) > maxTicks {
			step = int(math.Ceil(max / maxTicks))
		}
		ticks := make([]plot.Tick, 0, int(max)/step+2)
		for i := range int(max) + 1 {
			if i%step == 0 {
				ticks = append(ticks, plot.Tick{Value: float64(i), Label: fmt.Sprintf("%d", i)})
			} else {
				ticks = append(ticks, plot.Tick{Value: float64(i)})
			}
		}
		return ticks
	})
	p.Y.Min = 0
	p.Y.Max = 100
	pts := make(plotter.XYs, len(qstat)+1)
	pts[0].X = 0
	pts[0].Y = 100
	for i, qscore := range qstat {
		pts[i+1].X = float64(i + 1)
		pts[i+1].Y = 100 - qscore
	}
	line, points, err := plotter.NewLinePoints(pts)
	if err != nil {
		return err
	}
	line.Color = plotLineColor
	line.Dashes = []vg.Length{vg.Points(6), vg.Points(3)}
	points.Color = plotLineColor
	points.Shape = plotMarkerShap
	points.Radius = vg.Points(4)
	p.Add(line, points)
	return p.Save(plotW, plotH, fmt.Sprintf("%s.png", prefix))
}

// Write csv file containing reticulation branch scores to stdout
func WriteRetScoresToCSV(scores []*map[string]float64, names []string) error {
	branchNames := make([]string, 0)
	for k := range *scores[0] {
		branchNames = append(branchNames, k)
	}
	slices.SortFunc(branchNames, func(a, b string) int {
		if diff := len(a) - len(b); diff != 0 {
			return diff
		}
		return strings.Compare(a, b)
	})
	data := make([][]string, len(scores)+1)
	data[0] = append([]string{"gene"}, branchNames...)
	for i, row := range scores {
		data[i+1] = []string{names[i]}
		for _, br := range branchNames {
			data[i+1] = append(data[i+1], strconv.FormatFloat((*row)[br], 'f', -1, 64))
		}
	}
	writer := csv.NewWriter(os.Stdout)
	defer writer.Flush()
	if err := writer.WriteAll(data); err != nil {
		return fmt.Errorf("error writing csv file: %s", err)
	}
	return nil
}
