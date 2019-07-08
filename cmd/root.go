package cmd

import (
	"C"
	"bufio"
	"fmt"
	"io"
	"os"
	"runtime"
	"sync"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io/fasta"
	"github.com/evolbioinfo/gotree/io/utils"
	"github.com/evolbioinfo/gotree/tree"

	"github.com/spf13/cobra"
)

var treepath string
var alignpath string
var threads int

// rootCmd represents the base command when called without any subcommands
var rootCmd = &cobra.Command{
	Use:   "goliketree",
	Short: "Compute lk of a tree",
	Long:  ``,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var alignmentFile *os.File
		var al align.Alignment
		var treeChan <-chan tree.Trees
		var treeFile io.Closer
		var treeReader *bufio.Reader
		// Wait all Go routines
		var wg sync.WaitGroup

		if alignpath == "none" {
			err = fmt.Errorf("Alignment file is mandatory")
			return
		}

		if alignmentFile, err = os.Open(alignpath); err != nil {
			return
		}

		if al, err = fasta.NewParser(alignmentFile).Parse(); err != nil {
			return
		}
		patternWeightsInt := al.Compress()

		if treeFile, treeReader, err = utils.GetReader(treepath); err != nil {
			return
		}
		defer treeFile.Close()
		treeChan = utils.ReadMultiTrees(treeReader, utils.FORMAT_NEWICK)

		if threads > runtime.NumCPU() {
			threads = runtime.NumCPU()
		}
		runtime.GOMAXPROCS(threads)
		wg.Add(threads)

		for p := 0; p < threads; p++ {
			go func() {
				var err2 error
				var lk float64
				var trees tree.Trees
				defer wg.Done()
				for trees = range treeChan {
					if trees.Err != nil {
						err = trees.Err
						return
					}
					if err2 != nil {
						err = err2
						return
					}
					if lk, err2 = computelk(trees.Tree, al, patternWeightsInt); err != nil {
						err = err2
						return
					} else {
						fmt.Printf("Tree %d: lk=%v\n", trees.Id, lk)
					}
				}
			}()
		}
		wg.Wait()

		return
	},
}

// Execute adds all child commands to the root command and sets flags appropriately.
// This is called by main.main(). It only needs to happen once to the rootCmd.
func Execute() {
	if err := rootCmd.Execute(); err != nil {
		fmt.Println(err)
		os.Exit(1)
	}
}

func init() {
	rootCmd.PersistentFlags().StringVarP(&treepath, "tree", "t", "stdin", "Input tree file")
	rootCmd.PersistentFlags().IntVarP(&threads, "threads", "p", 1, "Number of threads")
	rootCmd.PersistentFlags().StringVarP(&alignpath, "align", "a", "none", "Input align file (mandatory)")
}
