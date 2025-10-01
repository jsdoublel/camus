package infer

import gr "github.com/jsdoublel/camus/internal/graphs"

// traceback for node v if there is not an edge (stored in DP.Traceback struct field)
type trace interface {
	traceback() []gr.Branch // returns all branches in subnetwork
}

// traceback if there isn't a cycle
type noCycleTrace struct {
	prevs [2]*trace // previous subproblems
}

func (tr *noCycleTrace) traceback() []gr.Branch {
	if tr.prevs[0] == nil {
		return []gr.Branch{}
	}
	return append((*tr.prevs[0]).traceback(), (*tr.prevs[1]).traceback()...)
}

// stores backtrace information along cycle
type cycleTraceNode struct {
	sib *trace          // sibling node trace
	p   *cycleTraceNode // parent node trace
}

func (tr *cycleTraceNode) traceUp() []gr.Branch {
	result := (*tr.sib).traceback()
	if tr.p != nil {
		result = append(result, tr.p.traceUp()...)
	}
	return result
}

// stores traceback info for node v in there is a cycle
type cycleTrace struct {
	pathW      *cycleTraceNode // beginning of linked-list w path towards v
	pathU      *cycleTraceNode // beginning of linked-list u path towards v
	wDownTrace *trace          // trace below w
	uDownTrace *trace          // trace below u
	branch     gr.Branch       // branch forming cycle
}

func (tr *cycleTrace) traceback() []gr.Branch {
	result := append((*tr.wDownTrace).traceback(), tr.branch)
	if tr.uDownTrace != nil {
		result = append(result, (*tr.uDownTrace).traceback()...)
	}
	if tr.pathU != nil {
		result = append(result, tr.pathU.traceUp()...)
	}
	if tr.pathW != nil {
		result = append(result, tr.pathW.traceUp()...)
	}
	return result
}
