#!/usr/bin/env	python3

"""
Filter multiple sequence alignment and transform it into a profile matrix

By Veronica Mixao
@INSA
"""

import argparse
import textwrap
import sys
from Bio import SeqIO
import pandas
import numpy as np


# functions	----------

def align2mx(align, log):
	""" transform an alignment into a matrix 
	input: fa
	output: pandas dataframe
	"""
	
	info = []
	col_names = ["#ID"]
	index_val = []
	
	for seq_record in SeqIO.parse(align, "fasta"):
		index_val.append(seq_record.id)
		line = [seq_record.id] + list(seq_record.seq)
		info.append(line)
		size = len(line) - 1
	
	for pos in range(1,size + 1):
		col_names.append("pos_" + str(pos))
	
	df = pandas.DataFrame(data = np.array(info), index = index_val, columns = col_names)
	
	return df


def clean_mx(align, mx, filters, log):
	""" remove samples according to metadata
	input: matrix
	output: filtered matrix
	"""
	
	mx = pandas.read_table(mx)
	sample_column = mx.columns[0]
	
	if "date" in mx.columns and "iso_week" not in mx.columns:
		index_no = mx.columns.get_loc("date")
		mx["date"] = pandas.to_datetime(mx["date"], errors = "coerce")
		year = mx["date"].dt.isocalendar().year
		week = mx["date"].dt.isocalendar().week
		mx["iso_year"] = year.astype(str)
		mx["iso_week"] = week.astype(str)
		mx["iso_date"] = year.astype(str).replace("<NA>", "-") + "W" + week.astype(str).replace("<NA>", "--")
		isoyear = mx.pop("iso_year")
		isoweek = mx.pop("iso_week")
		isodate = mx.pop("iso_date")
		mx.insert(index_no + 1, "iso_year", isoyear)
		mx.insert(index_no + 2, "iso_week", isoweek)
		mx.insert(index_no + 3, "iso_date", isodate)
		
	print("\tFiltering metadata for the following parameters: " + " & ".join(filters.split(";")))
	print("\tFiltering metadata for the following parameters: " + " & ".join(filters.split(";")), file = log)
	
	f = []
	if ";" in filters:
		for flt in filters.split(";"):
			f.append(flt)
	else:
		f.append(filters)

	for spec in f:
		col = spec.split(" ")[0]
		val = spec.split(" ")[1]
		cond = spec.split(" ")[2]
			
		if "," in cond:
			lst = cond.split(",")
			if val == "==":
				mx = mx[mx[col].isin(lst)]
			elif val == "!=":
				mx = mx[mx[col].isin(lst) == False]
		else:
			if col == "date":
				mx["date"] = mx["date"].astype("datetime64[ns]")
				if val == "==":
					mx = mx[mx["date"] == cond]  
				elif val == "!=":
					mx = mx[mx["date"] != cond] 
				elif val == ">":
					mx = mx[mx["date"] > cond] 
				elif val == ">=":
					mx = mx[mx["date"] >= cond] 
				elif val == "<=":
					mx = mx[mx["date"] <= cond] 
				elif val == "<":
					mx = mx[mx["date"] < cond]
						
			else:
				if val == "==":
					mx = mx[mx[col] == cond]
				elif val == "!=":
					mx = mx[mx[col] != cond]
				else:
					mx[col] = pandas.to_numeric(mx[col], errors='coerce')
					if val == ">":
						mx = mx[mx[col] > float(cond)]
					elif val == ">=":
						mx = mx[mx[col] >= float(cond)]
					elif val == "<":
						mx = mx[mx[col] < float(cond)]
					elif val == "<=":
						mx = mx[mx[col] >= float(cond)]
					mx[col] = mx[col].astype(int)
					mx[col] = mx[col].astype(str)

	samples = mx[sample_column].tolist()
	
	alignment_filtered = align[align[align.columns[0]].isin(samples)]
	
	return alignment_filtered
	

def clean_position(mx, log):
	""" remove gaps and conserved positions
	input: alignment
	output: alignment no gaps or conserved positions
	"""
	
	for col in mx.columns[1:]:
		values = mx[col].values.tolist()
		if "-" in values:
			mx = mx.drop(columns=col)
		elif len(set(values)) == 1:
			mx = mx.drop(columns=col)
		elif len(set(values)) == 2 and "N" in set(values):
			mx = mx.drop(columns=col)
	
	mx = mx.replace({"N": 0})
	
	return mx
	

# running the pipeline	----------

if __name__ == "__main__":
    
	# argument options
    
	parser = argparse.ArgumentParser(prog="partitioning_HC.py", formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent("""\
									###############################################################################             
									#                                                                             #
									#                          alignment_processing.py                            #
									#                                                                             #
									###############################################################################  
									                            
									alignment_processing.py cleans a multiple sequence alignment and transforms it
									into a profile matrix that can be used as input to partitioning_grapetree.py
									or to partitioning_HC.py.
									
									-----------------------------------------------------------------------------"""))
	
	group0 = parser.add_argument_group("Processing alignment", "Specifications to clean a multiple sequence alignment")
	group0.add_argument("-align", "--alignment", dest="alignment", default="", required=True, type=str, help="[MANDATORY] Input multiple sequence alignment (fasta format)")
	group0.add_argument("-o", "--output", dest="out", required=True, type=str, help="[MANDATORY] Tag for output file name")
	group0.add_argument("-m", "--metadata", dest="metadata", required=False, default="", type=str, help="[OPTIONAL] Metadata file in .tsv format to select the samples to reconstruct the minimum \
						spanning tree according to the '--filter' argument")
	group0.add_argument("-f", "--filter", dest="filter_column", required=False, default="", help="[OPTIONAL] Filter for metadata columns to select the samples of the allele matrix that must \
						be used for tree reconstruction. This must be specified within quotation marks in the following format 'column< >operation< >condition' (e.g. 'country == Portugal'). When \
						more than one condition is specified for a given column, they must be separated with commas (e.g 'country == Portugal,Spain,France'). When filters include more than one \
						column, they must be separated with semicolon (e.g. 'country == Portugal,Spain,France;date > 2018-01-01;date < 2022-01-01'). White spaces are important in this argument, \
						so, do not leave spaces before and after commas/semicolons.")

	args = parser.parse_args()


	# starting logs	----------

	log_name = args.out + ".log"
	log = open(log_name, "a+")
	
	print("\n-------------------- alignment_processing.py --------------------\n")
	print("\n-------------------- alignment_processing.py --------------------\n", file = log)
	print(" ".join(sys.argv))
	print(" ".join(sys.argv), file = log)
	
	
	# transform alignment into matrix	----------
	
	print("\nGetting the matrix of the alignment...\n")
	print("\nGetting the matrix of the alignment...\n", file = log)
	
	mx = align2mx(args.alignment, log)
	
	
	# remove samples according to metadata
	
	if args.metadata != "" and args.filter_column != "":
		print("Filtering the allele matrix...")
		print("Filtering the allele matrix...", file = log)
		mx = clean_mx(mx, args.metadata, args.filter_column, log)
	
	elif args.metadata != "" and args.filter_column == "":
		print("Metadata file was provided but no filter was found... I am confused :-(")
		print("Metadata file was provided but no filter was found... I am confused :-(", file = log)
		sys.exit()

	elif args.metadata == "" and args.filter_column != "":
		print("Metadata file was not provided but a filter was found... I am confused :-(")
		print("Metadata file was not provided but a filter was found... I am confused :-(", file = log)
		sys.exit()
	
	
	# remove gaps
	
	print("\nRemoving conserved positions and positions with gaps...\n")
	print("\nRemoving conserved positions and positions with gaps...\n", file = log)
	
	mx = clean_position(mx, log)
	
	mx.to_csv(args.out + "_profile.tsv", index = False, header=True, sep ="\t")

print("\nalignment_processing.py is done!")
print("\nalignment_processing.py is done!", file = log)

log.close()

	
	
