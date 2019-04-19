package cmd

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
	"sync"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/models/dna"
	"github.com/evolbioinfo/gotree/tree"
	"gonum.org/v1/gonum/mat"
)

var beagleOpNone = C.int(C.BEAGLE_OP_NONE)
var mux sync.Mutex // To initialize beagle instance

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

func convertToDoubleArr(input []int) *C.double {
	a := make([]C.double, len(input))
	for i, v := range input {
		a[i] = C.double(v)
	}
	return &a[0]
}

func convertFloatToDoubleArr(input []float64) *C.double {
	a := make([]C.double, len(input))
	for i, v := range input {
		a[i] = C.double(v)
	}
	return &a[0]
}

func convertDenseToDoubleArr(input *mat.Dense) *C.double {
	r, c := input.Dims()
	a := make([]C.double, r*c)
	n := 0
	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			a[n] = C.double(input.At(i, j))
			n++
		}
	}
	return &a[0]
}

func makeBeagleInstance(t *tree.Tree, alignment align.Alignment, patternWeightsInt []int) (instance C.int, err error) {
	mux.Lock()
	tipCount := len(t.Tips())
	partialsBufferCount := tipCount - 1
	compactBufferCount := tipCount
	stateCount := 4
	patternCount := alignment.Length()
	eigenBufferCount := 1
	matrixBufferCount := 2*tipCount - 1
	categoryCount := 1
	scaleBufferCount := 0
	resourceCount := 0
	var preferenceFlags C.long
	var requirementFlags C.long
	var returnInfo C.BeagleInstanceDetails
	var evec, ivec *mat.Dense
	var eval []float64

	instance = C.beagleCreateInstance(
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
		err = fmt.Errorf("Failed to obtain BEAGLE instance")
		return
	}

	patternWeights := convertToDoubleArr(patternWeightsInt)
	C.beagleSetPatternWeights(instance, patternWeights)

	freqs := filledDoubleArr(4, 0.25)
	C.beagleSetStateFrequencies(instance, 0, freqs)

	weights := filledDoubleArr(1, 1.)
	C.beagleSetCategoryWeights(instance, 0, weights)
	rates := filledDoubleArr(1, 1.)
	C.beagleSetCategoryRates(instance, rates)

	// evec := []C.double{1.0, 2.0, 0.0, 0.5, 1.0, -2.0, 0.5, 0.0, 1.0, 2.0, 0.0, -0.5, 1.0, -2.0, -0.5, 0.0}
	// ivec := []C.double{0.25, 0.25, 0.25, 0.25, 0.125, -0.125, 0.125, -0.125, 0.0, 1.0, 0.0, -1.0, 1.0, 0.0, -1.0, 0.0}
	// eval := []C.double{0.0, -1.3333333333333333, -1.3333333333333333, -1.3333333333333333}
	m := dna.NewGTRModel()
	m.InitModel(1., 1., 1., 1., 1., 1., 1./4., 1./4., 1./4., 1./4.)
	if eval, ivec, evec, err = m.Eigens(); err != nil {
		return
	}
	C.beagleSetEigenDecomposition(instance, 0,
		convertDenseToDoubleArr(evec),
		convertDenseToDoubleArr(ivec),
		convertFloatToDoubleArr(eval))
	mux.Unlock()
	return
}

func computelk(t *tree.Tree, alignment align.Alignment, patternWeightsInt []int) (loglk float64, err error) {
	stateMap := stateMapOfAlignment(alignment)

	var instance C.int
	if instance, err = makeBeagleInstance(t, alignment, patternWeightsInt); err != nil {
		return
	}
	defer C.beagleFinalizeInstance(instance)

	tipCount := len(t.Tips())
	edgeCount := 2*tipCount - 1

	if alignment.NbSequences() != tipCount {
		err = fmt.Errorf("Number of sequences doesn't match the number of tips.")
		return
	}

	edgeLengths := make([]C.double, edgeCount)
	nodeIndices := make([]C.int, edgeCount)

	nextID := 0
	// We go through the tips first because they need to have the lowest ids.
	t.PostOrder(func(cur *tree.Node, prev *tree.Node, parentEdge *tree.Edge) bool {
		if cur.Tip() {
			cur.SetId(nextID)
			nextID = nextID + 1
			nodeIndices[cur.Id()] = C.int(cur.Id())
			C.beagleSetTipStates(instance, C.int(cur.Id()), stateMap[cur.Name()])
			edgeLengths[cur.Id()] = C.double(parentEdge.Length())
		}
		return true
	})

	t.PostOrder(func(cur *tree.Node, prev *tree.Node, parentEdge *tree.Edge) bool {
		if !cur.Tip() {
			cur.SetId(nextID)
			nextID = nextID + 1
			nodeIndices[cur.Id()] = C.int(cur.Id())

			if parentEdge != nil {
				edgeLengths[cur.Id()] = C.double(parentEdge.Length())
			}
		}
		return true
	})

	C.beagleUpdateTransitionMatrices(instance,
		0, // eigenIndex
		(*C.int)(&nodeIndices[0]), // probabilityIndices
		nil, // firstDerivativeIndices
		nil, // secondDerivativeIndices
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

	t.PostOrder(func(cur *tree.Node, prev *tree.Node, parentEdge *tree.Edge) bool {
		if cur.Tip() == false {
			var leftChildID C.int
			var rightChildID C.int
			neigh := cur.Neigh()
			// For now, we assume a bifurcating root
			if cur == t.Root() {
				leftChildID = C.int(neigh[0].Id())
				rightChildID = C.int(neigh[1].Id())
			} else {
				if cur.Nneigh() != 3 {
					err = fmt.Errorf("Internal node doesn't have degree 3.")
					return false
				}
				if neigh[0] != prev {
					err = fmt.Errorf("Neighbors are not ordered as expected.")
					return false
				}
				leftChildID = C.int(neigh[1].Id())
				rightChildID = C.int(neigh[2].Id())
			}
			// fmt.Println("id:", cur.Id(), leftChildID, rightChildID)
			operations = append(operations, C.makeOperation(C.int(cur.Id()), beagleOpNone, beagleOpNone, leftChildID, leftChildID, rightChildID, rightChildID))
		}
		return true
	})
	if err != nil {
		return
	}
	if len(operations) == 0 {
		err = fmt.Errorf("No operations to do!")
		return
	}

	C.beagleUpdatePartials(instance,
		(*C.BeagleOperation)(&operations[0]),
		C.int(len(operations)),
		beagleOpNone)

	var logLp C.double
	rootIndex := [1]C.int{C.int(t.Root().Id())}
	categoryWeightIndex := [1]C.int{0}
	stateFrequencyIndex := [1]C.int{0}
	cumulativeScaleIndex := [1]C.int{beagleOpNone}

	C.beagleCalculateRootLogLikelihoods(instance,
		(*C.int)(&rootIndex[0]),
		(*C.int)(&categoryWeightIndex[0]),
		(*C.int)(&stateFrequencyIndex[0]),
		(*C.int)(&cumulativeScaleIndex[0]),
		1,
		&logLp)

	loglk = float64(logLp)

	return
}
