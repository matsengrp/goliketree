package cmd

import (
	"fmt"

	"github.com/spf13/cobra"
)

var Version string = "none"

// versionCmd represents the version command
var versionCmd = &cobra.Command{
	Use:   "version",
	Short: "Prints version",
	Run: func(cmd *cobra.Command, args []string) {
		fmt.Printf("goliketree %s\n",Version)
	},
}

func init() {
	rootCmd.AddCommand(versionCmd)
}
