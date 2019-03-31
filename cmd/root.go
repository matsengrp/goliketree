package cmd

import (
	"bufio"
	"fmt"
	"io"
	"os"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io/fasta"
	"github.com/evolbioinfo/gotree/io/utils"
	"github.com/evolbioinfo/gotree/tree"

	"github.com/spf13/cobra"
)

var treepath string
var alignpath string

func treesOfPath(treePath string) (treeSlice []*tree.Tree, err error) {
	var trees tree.Trees
	var treeChan <-chan tree.Trees
	var treeFile io.Closer
	var treeReader *bufio.Reader

	treeSlice = make([]*tree.Tree, 0, 10)

	if treeFile, treeReader, err = utils.GetReader(treePath); err != nil {
		panic(err)
	}
	defer treeFile.Close()
	treeChan = utils.ReadMultiTrees(treeReader, utils.FORMAT_NEWICK)
	for trees = range treeChan {
		if trees.Err != nil {
			err = trees.Err
			return
		}
		treeSlice = append(treeSlice, trees.Tree)
	}

	return
}

// rootCmd represents the base command when called without any subcommands
var rootCmd = &cobra.Command{
	Use:   "goliketree",
	Short: "Compute lk of a tree",
	Long:  ``,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var alignmentFile *os.File
		var trees []*tree.Tree
		var align align.Alignment

		if alignpath == "none" {
			err = fmt.Errorf("Alignment file is mandatory")
			return
		}

		if alignmentFile, err = os.Open(alignpath); err != nil {
			return
		}

		if align, err = fasta.NewParser(alignmentFile).Parse(); err != nil {
			return
		}

		if trees, err = treesOfPath(treepath); err != nil {
			return
		}
		err = computelk(trees, align)
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
	rootCmd.PersistentFlags().StringVarP(&alignpath, "align", "a", "none", "Input align file (mandatory)")
}
