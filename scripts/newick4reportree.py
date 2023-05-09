#!/usr/bin/env	python3

"""
Convert a Newick file into a newick format compatible with ReporTree

By Veronica Mixao
@INSA
"""

import argparse
import textwrap
import ete3

def read_tree(tree, tree_format):
	t = ete3.Tree(tree, format = tree_format)
	
	return t

def print_tree(tree, out):
	tree.write(outfile = str(out) + ".nwk")
	
def main():
    
	# argument options	----------
    
	parser = argparse.ArgumentParser(prog="newick4reportree.py", formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent("""\
									###############################################################################             
									#                                                                             #
									#                            newick4reportree.py                              #
									#                                                                             #
									############################################################################### 
									                            
									Convert a newick file into a Newick format compatible with ReporTree :-)
									                  
									-------------------------------------------------------------------------------"""))
	
	## parameters
	
	parser.add_argument("-t", "--tree", dest="tree", required=True, type=str, help="Input tree in Newick format")
	parser.add_argument("-o", "--out", dest="output", required=True, type=str, help="Tag for output file")
	parser.add_argument("-f", "--format", dest="tree_format", required=True, type=int, help="Format of the Newick file provided as input (check http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#reading-and-writing-newick-trees)")
	
	args = parser.parse_args()
	
	## pipeline
	
	t = read_tree(args.tree, args.tree_format)
	print_tree(t, args.output)
	
if __name__ == "__main__":
    main()
