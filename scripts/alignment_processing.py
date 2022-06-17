#!/usr/bin/env	python3

"""
Filter multiple sequence alignment and transform it into a profile matrix

By Veronica Mixao
@INSA
"""

import sys
import os
import argparse
import textwrap
from Bio import SeqIO, AlignIO
import pandas
import numpy as np


# functions	----------

def get_ref_coords(alignment, reference):
	""" get a dictionary with d[alignment position] = reference position
	input: alignment
	output: dictionary
	"""
	
	coords = {}
	i = 0 # ref
	j = 0 # del count
	k = 0 # align
	
	for record in alignment:
		if record.id == reference:
			seq = record.seq
	
	for nucl in seq:
		k += 1
		if nucl == "-":
			if j == 0:
				j = 1 # del count
			else:
				j += 1
			code = str(i) + "." + str(j)
			coords[k] = code
		else:
			i += 1
			j = 0
			coords[k] = i
			
	return coords
	
		
def filter_align(align, mx, filters, ref, log):
	""" remove samples according to metadata
	input: alignment
	output: filtered alignment
	"""
	
	mx = pandas.read_table(mx, dtype = str)
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
	
	flt_alignment = []
	for record in align:
		if record.id in samples and record.id != ref:
			flt_alignment.append(record)
		elif record.id not in samples and record.id == ref and ref != "none":
			flt_alignment.append(record)
	
	return flt_alignment
	

def core2mx(core, vcf, coords, log):
	""" transform the final alignment into a matrix
	input: alignment
	out: dataframe"""
	
	align_len = len(core[0].seq)
	
	# column names
	vcf_file = pandas.read_table(vcf)
	if coords != "":
		vcf_file["POS"] = vcf_file["POS"].replace(coords)
	
	# chunk details
	start = 0
	chunk_size = 100
	end = start + chunk_size
	
	# starting final dataframe
	sample_ids = []
	for record in core:
		sample_ids.append(record.id)
	final_df = pandas.DataFrame(data = sample_ids, index = sample_ids)		
	
	while end <= align_len + chunk_size:
		info = []
		col_names = []
		for record in core[:,start:end]:
			info.append(list(str(record.seq)))
		
		for pos in range(start,start + len(info[0])):
			col_names.append(str(pos + 1))
		
		df = pandas.DataFrame(data = np.array(info), index = sample_ids)
		final_df = pandas.concat([final_df, df], axis=1, join="inner")
		start += chunk_size
		end += chunk_size
	
	col_id = ["#ID"]
	pos_id = list(vcf_file["POS"])
	final_colnames = col_id + pos_id
	final_df.columns = final_colnames
	
	return final_df
	
		
def clean_mx(mx, gaps, log):
	""" clean alignment matrix
	input: matrix
	output: matrix
	"""
	
	align_len = len(mx.columns)
	
	df = mx.replace({"N": "0", "a": "A", "c": "C", "t": "T", "g": "G"})
	
	allowed_values = ["-", "0", "A", "T", "C", "G"]
	df[~df[df.columns[1:]].isin(allowed_values)] = "0"
		
	for col in df.columns:
		values = df[col].values.tolist()
		if len(set(values)) == 1:
			df = df.drop(columns=col)
		elif "-" in values and not gaps:
			df = df.drop(columns=col)
		elif len(set(values)) == 2 and "0" in set(values):
			df = df.drop(columns=col)
		elif len(set(values)) == 2 and "-" in set(values):
			mx = mx.drop(columns=col)
		elif len(set(values)) == 3 and "-" in set(values) and "0" in set(values):
			mx = mx.drop(columns=col)
	print("\tAlignment 1st round comprises " + str(len(df.index)) + " samples and " + str(len(df.columns) - 1) + " positions.")
	print("\tAlignment 1st round comprises " + str(len(df.index)) + " samples and " + str(len(df.columns) - 1) + " positions.", file = log)
	
	return df
	

def clean_position(mx, atcg, gaps, log):
	""" remove gaps and conserved positions
	input: alignment
	output: alignment no gaps or conserved positions
	"""
	
	mx = mx.replace({"N": "0", "a": "A", "c": "C", "t": "T", "g": "G"})
	
	allowed_values = ["-", "0", "A", "T", "C", "G"]
	mx[~mx[mx.columns[1:]].isin(allowed_values)] = "0"
	print("\tSites with < " + str(atcg*100) + "% ATCG will also be removed...")
	print("\tSites with < " + str(atcg*100) + "% ATCG will also be removed...", file = log)
	for col in mx.columns[1:]:
		values = mx[col].values.tolist()
		if "-" in values and not gaps:
			mx = mx.drop(columns=col)
		elif len(set(values)) == 1:
			mx = mx.drop(columns=col)
		elif len(set(values)) == 2 and "0" in set(values):
			mx = mx.drop(columns=col)
		elif (len(values)-values.count("0"))/len(values) < atcg:
			mx = mx.drop(columns=col)
		elif len(set(values)) == 2 and "-" in set(values):
			mx = mx.drop(columns=col)
		elif len(set(values)) == 3 and "-" in set(values) and "0" in set(values):
			mx = mx.drop(columns=col)
			
	print("\tCleaned alignment comprises " + str(len(mx.index)) + " samples and " + str(len(mx.columns) - 1) + " positions.")
	print("\tCleaned alignment comprises " + str(len(mx.index)) + " samples and " + str(len(mx.columns) - 1) + " positions.", file = log)

	return mx
	

def rm_ns(mx, ATCG_content, out, log):
	""" remove samples with <= ATCG_content from the alignment
	input = mx
	output = mx filtered """
	
	report_mx = {}
	len_align = len(mx.columns) - 1
	
	report_mx["samples"] = mx[mx.columns[0]]
	report_mx["missing"] = mx.isin(["0"]).sum(axis=1)
	report_mx["called"] = len_align - mx.isin(["0"]).sum(axis=1)	
	report_mx["pct_called"] = (len_align - mx.isin(["0"]).sum(axis=1)) / len_align

	report_df = pandas.DataFrame(report_mx)
	flt_report = report_df[report_df["pct_called"] > float(ATCG_content)]
	pass_samples = flt_report["samples"].values.tolist()
	
	removed_samples = set(mx[mx.columns[0]].values.tolist()) - set(pass_samples)
	
	print("\tFrom the " + str(len(mx[mx.columns[0]].values.tolist())) + " samples, " + str(len(pass_samples)) + " were kept in the alignment matrix.")
	print("\tFrom the " + str(len(mx[mx.columns[0]].values.tolist())) + " samples, " + str(len(pass_samples)) + " were kept in the alignment matrix.", file = log)
	
	if len(removed_samples) > 0:
		with open(out + "_samples_excluded.txt", "w+") as output:
			for sample in removed_samples:
				print(sample, file = output)
	
	report_df.to_csv(out + "_pos_report.tsv", index = False, header=True, sep ="\t")
	
	final_mx = mx[mx[mx.columns[0]].isin(pass_samples)]
	
	return final_mx


def df2fa(mx, output_file):
	""" convert a dataframe into a fasta file
	input: data frame
	output: fasta """
	
	mx = mx.replace({"0": "N"})
	seq_names = mx[mx.columns[0]].values.tolist()
	sequences = mx[mx.columns[1:]].values.tolist()
	
	with open(output_file, "w+") as out:
		for i in range(0,len(seq_names)):
			print(">" + str(seq_names[i]), file = out)
			print("".join(sequences[i]), file = out)


def rm_ref(align, ref):
	""" remove ref
	input: alignment
	output: alignment """
	
	flt_alignment = []
	for record in align:
		if record.id != ref:
			flt_alignment.append(record)
	
	return flt_alignment
		

def pos_int(mx, outfile):
	""" Generate a file with positions of interest
	input: dataframe
	output: single-column file """
	
	with open(outfile, "w+") as out:
		print("POS", file = out)
		for pos in mx.columns[1:]:
			print(pos, file = out)


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
									
									Conserved positions and those with 1 gaps will always be removed (unless you
									specify the argument '--keep-gaps').							
									
									Note: This script takes advantage of SNP-sites. Do not forget to cite its 
									authors as well!
									
									-----------------------------------------------------------------------------"""))
	
	group0 = parser.add_argument_group("Processing alignment", "Specifications to clean a multiple sequence alignment")
	group0.add_argument("-align", "--alignment", dest="alignment", default="", required=True, type=str, help="[MANDATORY] Input multiple sequence alignment (fasta format)")
	group0.add_argument("-o", "--output", dest="out", required=True, type=str, help="[MANDATORY] Tag for output file name")
	group0.add_argument("--sample-ATCG-content", dest="ATCG_content", required=False, default=0.0, help="[OPTIONAL: Useful to remove samples with excess of missing data] Minimum proportion (\
						between 0 and 1) of ATCG in informative sites of the alignment per sample (e.g. '--ATCG-content 1.0' will only keep samples without N's or any non-ATCG code in \
						informative sites) NOTE: This argument works on samples (i.e. rows) and will be applied after gap and conserved positions removal and before '--site-ATCG-content'! \
						[default: 0 - keep all samples]")
	group0.add_argument("--site-ATCG-content", dest="N_content", required=False, default=0.0, help="[OPTIONAL: Useful to remove informative sites with excess of missing data] Minimum proportion \
						(between 0 and 1) of ATCG per informative site of the alignment (e.g. '--site-ATCG-content 1.0' will only keep positions without N's or any non-ATCG code, i.e. a core \
						SNP alignment) NOTE: This argument works on alignment positions (i.e. columns)! [default: 0.0 - content of missing data is not considered during alignment cleaning]")
	group0.add_argument("--remove-reference", dest="remove_ref", required=False, action="store_true", help="Set only if you want to remove the reference sequence of the alignment (reference \
						name must be provided with the argument '--reference').")
	group0.add_argument("--keep-gaps", dest="keep_gaps", required=False, action="store_true", help="Set only if you want that columns with gaps are not removed just because they have a gap.")	
	group0.add_argument("--use-alignment-coords", dest="use_align", required=False, action="store_true", help="Set only if you want that column names in the final matrix represent the initial \
						alignment coordinates. Note: Depending on the alignment size, this argument can make alignment processing very slow!")	
	group0.add_argument("--use-reference-coords", dest="use_ref", required=False, action="store_true", help="Set only if you want that column names in the final matrix represent the reference \
						coordinates (reference name must be provided with the argument '--reference') Note: Depending on the alignment size, this argument can make alignment processing very \
						slow!")	
	group0.add_argument("-r", "--reference", dest="reference", required=False, type=str, default = "none", help="[OPTIONAL] Name of reference sequence. Required if '--remove-reference' and/or \
						'--use-reference-coords' specified.")
	group0.add_argument("-m", "--metadata", dest="metadata", required=False, default="", type=str, help="[OPTIONAL] Metadata file in .tsv format to apply sample subset.")
	group0.add_argument("-f", "--filter", dest="filter_column", required=False, default="", help="[OPTIONAL] Filter for metadata columns to select the samples of the alignment that must \
						be included in the matrix. This must be specified within quotation marks in the following format 'column< >operation< >condition' (e.g. 'country == Portugal'). When \
						more than one condition is specified for a given column, they must be separated with commas (e.g 'country == Portugal,Spain,France'). When filters include more than one \
						column, they must be separated with semicolon (e.g. 'country == Portugal,Spain,France;date > 2018-01-01;date < 2022-01-01'). White spaces are important in this argument, \
						so, do not leave spaces before and after commas/semicolons.")
	group0.add_argument("--get-positions-interest", dest="pos_int", required=False, action="store_true", help="Set only if you want a file with the final positions of interest.")
	
	args = parser.parse_args()


	# starting logs	----------
	
	log_name = args.out + ".log"
	log = open(log_name, "a+")
	
	print("\n-------------------- alignment_processing.py --------------------\n")
	print("\n-------------------- alignment_processing.py --------------------\n", file = log)
	print(" ".join(sys.argv))
	print(" ".join(sys.argv), file = log)
	
	if args.remove_ref and args.reference == "none":
		print("'--remove-reference' specified, but I do not know what is the reference id :-(")
		print("'--remove-reference' specified, but I do not know what is the reference id :-(", file = log)
		sys.exit()
	elif args.use_ref and args.reference == "none":
		print("'--use-reference-coords' specified, but I do not know what is the reference id :-(")
		print("'--use-reference-coords' specified, but I do not know what is the reference id :-(", file = log)
		sys.exit()
		
		
	# read alignment
	
	print("Loading the alignment...")
	print("Loading the alignment...", file = log)
	
	alignment = AlignIO.read(args.alignment, "fasta")
	align_len = len(alignment[0].seq)
	print("\tLoaded " + str(len(alignment)) + " samples.")
	print("\tLoaded " + str(len(alignment)) + " samples.", file = log)
	print("\tInitial alignment length: " + str(align_len))
	print("\tInitial alignment length: " + str(align_len), file = log)

	# getting coordinate correspondence
	if args.reference != "none" and args.use_ref:
		coords = get_ref_coords(alignment, args.reference)
		
	# remove samples according to metadata
	
	if args.metadata != "" and args.filter_column != "":
		print("Filtering the alignment...")
		print("Filtering the alignment...", file = log)
		alignment = filter_align(alignment, args.metadata, args.filter_column, args.reference, log)
		SeqIO.write(alignment, "tmp.fasta", "fasta")
		alignment = AlignIO.read("tmp.fasta", "fasta")
		os.system("rm tmp.fasta")
	elif args.metadata != "" and args.filter_column == "":
		print("Metadata file was provided but no filter was found... I am confused :-(")
		print("Metadata file was provided but no filter was found... I am confused :-(", file = log)
		sys.exit()

	elif args.metadata == "" and args.filter_column != "":
		print("Metadata file was not provided but a filter was found... I am confused :-(")
		print("Metadata file was not provided but a filter was found... I am confused :-(", file = log)
		sys.exit()
	
	# run snp-sites to make the alignment shorter (only if '--use-reference-coords' and '--use-alignment-coords' are not set
	
	if float(args.ATCG_content) == 0.0 and float(args.N_content) == 1.0: # do not care about samples N content and keep only sites with ATCG
		if args.remove_ref:
			print("Removing reference sequence: " + str(args.reference))
			print("Removing reference sequence: " + str(args.reference), file = log)
			alignment = rm_ref(alignment, args.reference)
		print("Trimming the alignment with SNP-sites:")
		print("Trimming the alignment with SNP-sites:", file = log)
		print("\tsnp-sites -c tmp.fasta > tmp_flt.fasta")
		print("\tsnp-sites -c tmp.fasta > tmp_flt.fasta", file = log)
		SeqIO.write(alignment, "tmp.fasta", "fasta")
		os.system("snp-sites -c tmp.fasta > tmp_flt.fasta")
		os.system("snp-sites -v -c tmp.fasta > tmp_flt.vcf")
		os.system("rm tmp.fasta")
		alignment = AlignIO.read("tmp_flt.fasta", "fasta")
		os.system("grep -v '##' tmp_flt.vcf >  tmp.vcf")
		os.system("rm tmp_flt.fasta tmp_flt.vcf")
		print("\tAlignment length after SNP-sites: " + str(len(alignment[0].seq)))
		print("\tAlignment length after SNP-sites: " + str(len(alignment[0].seq)), file = log)
		
		# alignment matrix
		
		print("Getting the alignment matrix...")
		print("Getting the alignment matrix...", file = log)
		if args.use_ref:
			mx = core2mx(alignment, "tmp.vcf", coords, log)
			os.system("rm tmp.vcf")
		else:
			mx = core2mx(alignment, "tmp.vcf", "", log)
			os.system("rm tmp.vcf")
	
	else: # run snp-sites with option -c is not viable
		if args.remove_ref:
			print("Removing reference sequence: " + str(args.reference))
			print("Removing reference sequence: " + str(args.reference), file = log)
			alignment = rm_ref(alignment, args.reference)
		print("Trimming the alignment with SNP-sites:")
		print("Trimming the alignment with SNP-sites:", file = log)
		print("\tsnp-sites tmp.fasta > tmp_flt.fasta")
		print("\tsnp-sites tmp.fasta > tmp_flt.fasta", file = log)
		SeqIO.write(alignment, "tmp.fasta", "fasta")
		os.system("snp-sites tmp.fasta > tmp_flt.fasta")
		os.system("snp-sites -v tmp.fasta > tmp_flt.vcf")
		os.system("rm tmp.fasta")
		alignment = AlignIO.read("tmp_flt.fasta", "fasta")
		os.system("grep -v '##' tmp_flt.vcf >  tmp.vcf")
		os.system("rm tmp_flt.fasta tmp_flt.vcf")
		print("\tAlignment length after SNP-sites: " + str(len(alignment[0].seq)))
		print("\tAlignment length after SNP-sites: " + str(len(alignment[0].seq)), file = log)
		
		# alignment matrix
		
		print("Getting the alignment matrix...")
		print("Getting the alignment matrix...", file = log)
		if args.use_ref:
			mx = core2mx(alignment, "tmp.vcf", coords, log)
			os.system("rm tmp.vcf")
		else:
			mx = core2mx(alignment, "tmp.vcf", "", log)
			os.system("rm tmp.vcf")
		
		# additional cleaning
			
		print("Removing conserved positions and positions with gaps (1st round)...")
		print("Removing conserved positions and positions with gaps (1st round)...", file = log)
		mx = clean_mx(mx, args.keep_gaps, log)
				
		# assess number of ATCG's per sample
		
		print("Removing samples with <= " + str(float(args.ATCG_content)*100) + "% of ATCG in informative sites...")
		print("Removing samples with <= " + str(float(args.ATCG_content)*100) + "% of ATCG in informative sites...", file = log)
		mx_Ns = rm_ns(mx, args.ATCG_content, args.out, log)
		if len(mx_Ns.index) <= 1:
			print("Cannot proceed because " + str(len(mx_Ns.index)) + " samples were kept in the alignment!")
			print("Cannot proceed because " + str(len(mx_Ns.index)) + " samples were kept in the alignment!", file = log)
			sys.exit()

		print("Removing conserved positions and positions with gaps (2nd round)...")
		print("Removing conserved positions and positions with gaps (2nd round)...", file = log)
		mx = clean_position(mx_Ns, float(args.N_content), args.keep_gaps, log)
			
		if len(mx.index) <= 1:
			print("Cannot proceed because " + str(len(mx.index)) + " samples were kept in the alignment!")
			print("Cannot proceed because " + str(len(mx.index)) + " samples were kept in the alignment!", file = log)
			sys.exit()
		if len(mx.columns) <= 2:
			print("Cannot proceed because " + str(len(mx.columns)-1) + " sites were kept in the alignment!")
			print("Cannot proceed because " + str(len(mx.columns)-1) + " sites were kept in the alignment!", file = log)
			sys.exit()	
	
	# outputs
	
	mx.to_csv(args.out + "_align_profile.tsv", index = False, header=True, sep ="\t")
	df2fa(mx, args.out + "_align_profile.fasta")
	if args.pos_int:
		pos_int(mx, args.out + "_positions_of_interest.tsv")

print("\nalignment_processing.py is done!")
print("\nalignment_processing.py is done!", file = log)

log.close()
