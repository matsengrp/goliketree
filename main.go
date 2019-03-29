package main

// #cgo CFLAGS: -I/usr/local/include/libhmsbeagle-1
// #cgo LDFLAGS: -lhmsbeagle -L/usr/local/lib
//
// #include <libhmsbeagle/beagle.h>
//
// BeagleOperation makeOperation(
//     int destinationPartials,
//     int destinationScaleWrite,
//     int destinationScaleRead,
//     int child1Partials,
//     int child1TransitionMatrix,
//     int child2Partials,
//     int child2TransitionMatrix) {
//
// 	BeagleOperation op = {
// 		destinationPartials = destinationPartials,
// 		destinationScaleWrite = destinationScaleWrite,
// 		destinationScaleRead = destinationScaleRead,
// 		child1Partials = child1Partials,
// 		child1TransitionMatrix = child1TransitionMatrix,
// 		child2Partials = child2Partials,
// 		child2TransitionMatrix = child2TransitionMatrix,
// 	};
//
// 	return op;
// }
import "C"
import (
	"fmt"
	"log"
	"os"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io/fasta"
	"github.com/evolbioinfo/gotree/io/newick"
	"github.com/evolbioinfo/gotree/tree"
)

var BEAGLE_OP_NONE = C.int(C.BEAGLE_OP_NONE)

func getTable() [128]C.int {
	var table [128]C.int
	table[int('A')] = 0
	table[int('C')] = 1
	table[int('G')] = 2
	table[int('T')] = 3
	table[int('a')] = 0
	table[int('c')] = 1
	table[int('g')] = 2
	table[int('t')] = 3
	table[int('-')] = 4
	return table
}

func createStates(s []rune, table [128]C.int) *C.int {
	a := make([]C.int, len(s))
	for i := range s {
		a[i] = table[int(s[i])]
	}
	return (*C.int)(&a[0])
}

func stateMapOfAlignment(alignment align.Alignment) map[string]*C.int {
	m := make(map[string]*C.int)
	table := getTable()

	alignment.IterateChar(func(name string, sequence []rune) {
		m[name] = createStates(sequence, table)
	})
	return m
}

func filledDoubleArr(length int, value C.double) *C.double {
	a := make([]C.double, length)
	for i := range a {
		a[i] = value
	}
	return (*C.double)(&a[0])
}

func printNode(cur *tree.Node, prev *tree.Node) {
	fmt.Print(cur.Id(), ",")
}

func main() {
	alignmentPath := "hello.fasta"
	treePath := "hello.nwk"

	var alignmentFile *os.File
	alignmentFile, err := os.Open(alignmentPath)
	if err != nil {
		panic(err)
	}
	alignment, err := fasta.NewParser(alignmentFile).Parse()
	if err != nil {
		panic(err)
	}
	stateMap := stateMapOfAlignment(alignment)

	var t *tree.Tree
	var treeFile *os.File
	if treeFile, err = os.Open(treePath); err != nil {
		panic(err)
	}
	t, err = newick.NewParser(treeFile).Parse()
	if err != nil {
		panic(err)
	}

	tipCount := len(t.Tips())
	edgeCount := 2*tipCount - 1

	partialsBufferCount := tipCount - 1
	compactBufferCount := tipCount
	stateCount := 4
	patternCount := alignment.Length()
	eigenBufferCount := 1
	matrixBufferCount := edgeCount
	categoryCount := 1
	scaleBufferCount := 0
	resourceCount := 0
	var preferenceFlags C.long = 0
	var requirementFlags C.long = 0
	var returnInfo C.BeagleInstanceDetails

	instance := C.beagleCreateInstance(
		C.int(tipCount),
		C.int(partialsBufferCount),
		C.int(compactBufferCount),
		C.int(stateCount),
		C.int(patternCount),
		C.int(eigenBufferCount),
		C.int(matrixBufferCount),
		C.int(categoryCount),
		C.int(scaleBufferCount),
		nil, // resourceList,
		C.int(resourceCount),
		preferenceFlags,
		requirementFlags,
		&returnInfo)

	if instance < 0 {
		log.Fatal("Failed to obtain BEAGLE instance")
	}

	edgeLengths := make([]C.double, edgeCount)
	nodeIndices := make([]C.int, edgeCount)

	nextId := 0

	// We go through the tips first because they need to have the lowest ids.
	t.PostOrder(func(cur *tree.Node, prev *tree.Node, parentEdge *tree.Edge) {
		if cur.Tip() {
			cur.SetId(nextId)
			nextId = nextId + 1
			nodeIndices[cur.Id()] = C.int(cur.Id())
			C.beagleSetTipStates(instance, C.int(cur.Id()), stateMap[cur.Name()])
			edgeLengths[cur.Id()] = C.double(parentEdge.Length())
		}
	})

	t.PostOrder(func(cur *tree.Node, prev *tree.Node, parentEdge *tree.Edge) {
		if !cur.Tip() {
			cur.SetId(nextId)
			nextId = nextId + 1
			nodeIndices[cur.Id()] = C.int(cur.Id())

			if parentEdge != nil {
				edgeLengths[cur.Id()] = C.double(parentEdge.Length())
			}
		}
	})

	patternWeights := filledDoubleArr(patternCount, 1.)
	C.beagleSetPatternWeights(instance, patternWeights)

	freqs := filledDoubleArr(4, 0.25)
	C.beagleSetStateFrequencies(instance, 0, freqs)

	weights := filledDoubleArr(1, 1.)
	C.beagleSetCategoryWeights(instance, 0, weights)
	rates := filledDoubleArr(1, 1.)
	C.beagleSetCategoryRates(instance, rates)

	evec := []C.double{1.0, 2.0, 0.0, 0.5, 1.0, -2.0, 0.5, 0.0, 1.0, 2.0, 0.0, -0.5, 1.0, -2.0, -0.5, 0.0}
	ivec := []C.double{0.25, 0.25, 0.25, 0.25, 0.125, -0.125, 0.125, -0.125, 0.0, 1.0, 0.0, -1.0, 1.0, 0.0, -1.0, 0.0}
	eval := []C.double{0.0, -1.3333333333333333, -1.3333333333333333, -1.3333333333333333}
	C.beagleSetEigenDecomposition(instance, 0, &evec[0], &ivec[0], &eval[0])

	C.beagleUpdateTransitionMatrices(instance,
		0,                            // eigenIndex
		(*C.int)(&nodeIndices[0]),    // probabilityIndices
		nil,                          // firstDerivativeIndices
		nil,                          // secondDerivativeIndices
		(*C.double)(&edgeLengths[0]), // edgeLengths
		C.int(edgeCount))             // count

	operationCount := tipCount - 1
	operations := make([]C.BeagleOperation, 0, operationCount)

	// fmt.Println(nodeIndices)
	// fmt.Println(edgeLengths)
	// m := make([]C.double, 16)
	// for i := 0; i < edgeCount; i++ {
	// 	C.beagleGetTransitionMatrix(instance, C.int(i), (*C.double)(&m[0]))
	// 	fmt.Println(m)
	// }

	t.PostOrder(func(cur *tree.Node, prev *tree.Node, parentEdge *tree.Edge) {
		if cur.Tip() == false {
			var left_child_id C.int
			var right_child_id C.int
			neigh := cur.Neigh()
			// For now, we assume a bifurcating root
			if cur == t.Root() {
				left_child_id = C.int(neigh[0].Id())
				right_child_id = C.int(neigh[1].Id())
			} else {
				if cur.Nneigh() != 3 {
					panic("Internal node doesn't have degree 3.")
				}
				if neigh[0] != prev {
					panic("Neighbors are not ordered as expected.")
				}
				left_child_id = C.int(neigh[1].Id())
				right_child_id = C.int(neigh[2].Id())
			}
			// fmt.Println("id:", cur.Id(), left_child_id, right_child_id)
			operations = append(operations, C.makeOperation(C.int(cur.Id()), BEAGLE_OP_NONE, BEAGLE_OP_NONE, left_child_id, left_child_id, right_child_id, right_child_id))
		}
	})

	if len(operations) == 0 {
		panic("No operations to do!")
	}

	C.beagleUpdatePartials(instance,
		(*C.BeagleOperation)(&operations[0]),
		C.int(len(operations)),
		BEAGLE_OP_NONE)

	var logLp C.double
	rootIndex := [1]C.int{C.int(t.Root().Id())}
	categoryWeightIndex := [1]C.int{0}
	stateFrequencyIndex := [1]C.int{0}
	cumulativeScaleIndex := [1]C.int{BEAGLE_OP_NONE}

	C.beagleCalculateRootLogLikelihoods(instance,
		(*C.int)(&rootIndex[0]),
		(*C.int)(&categoryWeightIndex[0]),
		(*C.int)(&stateFrequencyIndex[0]),
		(*C.int)(&cumulativeScaleIndex[0]),
		1,
		&logLp)

	fmt.Println("logL =", logLp)
	fmt.Println("Woof!")
}
