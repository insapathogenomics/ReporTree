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
import pandas
import numpy as np
from datetime import date
import datetime as datetime
from Bio import SeqIO, AlignIO, Align, Alphabet

version = "1.4.0"
last_updated = "2023-12-11"

# functions	----------

def get_ref_coords(alignment, reference, sequence_corr, tag, pos_list):
	""" get a dictionary with d[alignment position] = reference position and a correspondence table if requested by the user
	input: alignment
	output: dictionary
	"""
	
	coords = {} #coords[alignment position] = reference position 
	coords_2_report = {} #coords_2_report[sequence][position] = alignment position 
	i = 0 # seq
	j = 0 # del count
	k = 0 # align
	
	corr_df = pandas.DataFrame()
	#corr = {} #corr[sample] = list(positions)
	corr_df["alignment_position"] = list(range(1, len(alignment[0].seq) + 1))
	i = 0
	for record in alignment:
		i += 1
		if record.id == reference:
			seq = record.seq
			name = "REFERENCE_" + record.id
			corr_info = []
			#corr["REFERENCE_" + record.id] = []
			coords_2_report["REFERENCE_" + record.id] = {}
	
			for nucl in seq:
				k += 1
				if nucl == "-":
					j += 1
					code = str(i) + "." + str(j)
					coords[k] = code
					#corr["REFERENCE_" + record.id].append(code)
					corr_info.append(code)
				else:
					i += 1
					j = 0
					coords[k] = i
					#corr["REFERENCE_" + record.id].append(str(i)  + " (" + nucl + ")")
					corr_info.append(str(i)  + " (" + nucl + ")")
					coords_2_report["REFERENCE_" + record.id][i] = k
		else:
			name = record.id
			corr_info = []
			if sequence_corr == "all":
				#corr[record.id] = []
				coords_2_report[record.id] = {}
				seq = record.seq
				counter = 0
				counter_align = 0
				gap_counter = 0
				for nucl in seq:
					counter_align += 1
					if nucl == "-":
						gap_counter += 1
						code = str(counter) + "." + str(gap_counter)
						#corr[record.id].append(code)
						corr_info.append(code)
					else:
						counter += 1
						#corr[record.id].append(str(counter) + " (" + nucl + ")")
						corr_info.append(str(counter)  + " (" + nucl + ")")
						coords_2_report[record.id][counter] = counter_align
						gap_counter = 0					
			elif "," in sequence_corr:
				if record.id in sequence_corr.split(","):
					#corr[record.id] = []
					coords_2_report[record.id] = {}
					seq = record.seq
					counter = 0
					counter_align = 0
					gap_counter = 0
					for nucl in seq:
						counter_align += 1
						if nucl == "-":
							gap_counter += 1
							code = str(counter) + "." + str(gap_counter)
							#corr[record.id].append(code)
							corr_info.append(code)
						else:
							counter += 1
							#corr[record.id].append(str(counter) + " (" + nucl + ")")
							corr_info.append(str(counter)  + " (" + nucl + ")")
							coords_2_report[record.id][counter] = counter_align
							gap_counter = 0
		corr_df[name] = corr_info	

	#corr_df = pandas.DataFrame(data = corr)

	if sequence_corr != "none":
		if pos_list != "none":
			pos_df = pandas.read_table(pos_list)
			pos_2_keep = []
			for sequence in pos_df.columns:
				for pos in pos_df[sequence]:
					if pos in coords_2_report[sequence].keys() and coords_2_report[sequence][pos]-1 not in pos_2_keep:
						pos_2_keep.append(coords_2_report[sequence][pos]-1)
			corr_df = corr_df.filter(items = pos_2_keep, axis=0)
		corr_df.to_csv(tag + "_align_position_correspondence.tsv", index = False, header=True, sep ="\t")
	
	return coords, corr_df
	
		
def filter_align(align, mx, filters, log):
	""" remove samples according to metadata
	input: alignment
	output: filtered alignment
	"""
	
	mx = pandas.read_table(mx, dtype = str)
	columns_names = [col.replace(" ", "_") for col in mx.columns]
	mx.replace(" ", "_", regex=True, inplace=True)
	mx.columns = columns_names
	sample_column = mx.columns[0]

	if "date" in mx.columns:
		index_no = mx.columns.get_loc("date")
		mx["date_original"] = mx["date"]
		date_original = mx.pop("date_original")
		mx.insert(index_no, "date_original", date_original)
		index_no = mx.columns.get_loc("date")
		if "year" in mx.columns:
			mx["year_original"] = mx["year"]
			year_original = mx.pop("year_original")
			mx.insert(index_no, "year_original", year_original)
			index_no = mx.columns.get_loc("date")
		mx["date"] = pandas.to_datetime(mx["date"], errors = "coerce")
		mx["year"] = mx["date"].dt.year
		year = mx.pop("year")
		mx.insert(index_no + 1, "year", year)
		index_no = mx.columns.get_loc("date")
		if "iso_week_nr" not in mx.columns and "iso_year" not in mx.columns and "iso_week" not in mx.columns:
			isoyear = mx["date"].dt.isocalendar().year
			isoweek = mx["date"].dt.isocalendar().week
			mx["iso_year"] = isoyear.astype(str)
			mx["iso_week_nr"] = isoweek.astype(str)
			mx["iso_week"] = isoyear.astype(str).replace("<NA>", "-") + "W" + isoweek.astype(str).replace("<NA>", "--").apply(lambda x: x.zfill(2))
			isoyear = mx.pop("iso_year")
			isoweek = mx.pop("iso_week_nr")
			isodate = mx.pop("iso_week")
			mx.insert(index_no + 2, "iso_year", isoyear)
			mx.insert(index_no + 3, "iso_week_nr", isoweek)
			mx.insert(index_no + 4, "iso_week", isodate)

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
		cond = " ".join(spec.split(" ")[2:])
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
			elif col == "iso_week":
				if "date" in mx.columns:
					year = cond.split("W")[0]
					week = cond.split("W")[1]
					cond = pandas.to_datetime(date.fromisocalendar(int(year), int(week), 1))
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
					print("\tCannot apply the 'iso_week' filter because column 'date' was not found in the metadata!")
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
		if record.id in samples:
			flt_alignment.append(record)

	print("\tFrom the " + str(len(align)) + " samples, " + str(len(flt_alignment)) + " were kept in the matrix...")
	print("\tFrom the " + str(len(align)) + " samples, " + str(len(flt_alignment)) + " were kept in the matrix...", file = log)

	return flt_alignment
	

def core2mx(core, vcf, coords, log):
	""" transform the final alignment into a matrix
	input: alignment
	out: dataframe 
	"""
	
	align_len = len(core[0].seq)
	
	# column names
	if vcf != "":
		vcf_file = pandas.read_table(vcf)
		if coords != "":
			vcf_file["POS"] = vcf_file["POS"].replace(coords)
		pos_id = list(vcf_file["POS"])
	else:
		pos_id = []
		for i in range(1,align_len + 1):
			if coords != "":
				pos_id.append(coords[i])
			else:
				pos_id.append(i)
		
	# chunk details
	start = 0
	chunk_size = 1000
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
	final_colnames = col_id + pos_id
	final_df.columns = final_colnames

	return final_df
	
		
def clean_mx(mx, gaps, all_gaps, missing_code, log):
    """ clean alignment matrix
    input: matrix
    output: matrix
    """
    
    align_len = len(mx.columns)
    df = mx.apply(lambda x: x.astype(str).str.upper())
    allowed_values = list(Alphabet.Gapped(Alphabet.IUPAC.ExtendedIUPACProtein, "-").letters)
    df = df.replace({missing_code: "0"})
    df[df.columns[0]] = mx[mx.columns[0]]

    df[~df[df.columns[1:]].isin(allowed_values)] = "0"
    
    final_df = pandas.DataFrame()
    start = 1
    chunk_size = 1000
    end = start + chunk_size
    final_df[df.columns[0]] = df[df.columns[0]]

    while end <= align_len + chunk_size:
        for col in df.columns[start:end]:
            append = True
            if len(df[col].unique()) == 1:
                append = False
            elif not gaps and not all_gaps and df[col].str.contains("-").sum() > 0:
                append = False
            elif len(df[col].unique()) == 2 and df[col].str.contains("0").sum() > 0:
                append = False
            elif len(df[col].unique()) == 2 and df[col].str.contains("-").sum() > 0 and not all_gaps:
                append = False
            elif len(df[col].unique()) == 3 and df[col].str.contains("-").sum() > 0 and df[col].str.contains("0").sum() > 0 and not all_gaps:
                append = False
            if append:
                final_df[col] = df[col]
        start += chunk_size
        end += chunk_size
			
		
    print("\tAlignment 1st round comprises " + str(len(df.index)) + " samples and " + str(len(final_df.columns) - 1) + " positions.")
    print("\tAlignment 1st round comprises " + str(len(df.index)) + " samples and " + str(len(final_df.columns) - 1) + " positions.", file = log)

    return final_df


def clean_position(mx, atcg, gaps, all_gaps, missing_code, log):
	""" remove gaps and conserved positions
	input: alignment
	output: alignment no gaps or conserved positions
	"""
	
	mx_id = mx[mx.columns[0]]
	mx = mx.apply(lambda x: x.astype(str).str.upper())
    
	allowed_values = list(Alphabet.Gapped(Alphabet.IUPAC.ExtendedIUPACProtein, "-").letters)
	mx = mx.replace({missing_code: "0"})
	mx[mx.columns[0]] = mx_id
	mx[~mx[mx.columns[1:]].isin(allowed_values)] = "0"
	
	print("\tSites with < " + str(atcg*100) + "% ATCG will also be removed...")
	print("\tSites with < " + str(atcg*100) + "% ATCG will also be removed...", file = log)
	
	final_df = pandas.DataFrame()
	start = 1
	chunk_size = 1000
	end = start + chunk_size
	final_df[mx.columns[0]] = mx[mx.columns[0]]

	while end <= align_len + chunk_size:
		for col in mx.columns[start:end]:
			append = True
			if mx[col].str.contains("-").sum() > 0 and not gaps and not all_gaps:
				append = False
			elif len(mx[col].unique()) == 1:
				append = False
			elif len(mx[col].unique()) == 2 and mx[col].str.contains("0").sum() > 0:
				append = False
			elif len(mx[col].unique()) == 2 and mx[col].str.contains("-").sum() > 0 and not all_gaps:
				append = False
			elif len(mx[col].unique()) == 3 and mx[col].str.contains("-").sum() > 0 and mx[col].str.contains("0").sum() > 0 and not all_gaps:
				append = False
			else:
				values = mx[col].values.tolist()
				if (len(values)-values.count("0"))/len(values) < atcg:
					append = False
			if append:
				final_df[col] = mx[col]
		start += chunk_size
		end += chunk_size
	
	print("\tCleaned alignment comprises " + str(len(final_df.index)) + " samples and " + str(len(final_df.columns) - 1) + " positions.")
	print("\tCleaned alignment comprises " + str(len(final_df.index)) + " samples and " + str(len(final_df.columns) - 1) + " positions.", file = log)

	return final_df
	

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
	if float(ATCG_content) != 1.0:
		flt_report = report_df[report_df["pct_called"] > float(ATCG_content)]
	else:
		flt_report = report_df[report_df["pct_called"] == float(ATCG_content)]
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
    
	parser = argparse.ArgumentParser(prog="alignment_processing.py", formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent("""\
									###############################################################################             
									#                                                                             #
									#                          alignment_processing.py                            #
									#                                                                             #
									###############################################################################  
									                            
									alignment_processing.py cleans a multiple sequence alignment and transforms it
									into a profile matrix that can be used as input to partitioning_grapetree.py
									or to partitioning_HC.py. It works with the IUPAC nomenclature!

									Gaps should always be indicated with "-".
									By default, "N" is interpreted as missing data. However, other codes can be 
									indicated with the '--missing-code' option.
									
									Conserved positions and those with 1 gaps will always be removed (unless you
									specify the argument '--keep-gaps' or '--keep-all-gaps').							
									
									Please avoid blank spaces in sample names, because only the first part of the
									name will be used!

									Note: This script takes advantage of SNP-sites. Do not forget to cite its 
									authors as well!
									
									-----------------------------------------------------------------------------"""))
	
	group0 = parser.add_argument_group("Processing alignment", "Specifications to clean a multiple sequence alignment (nucleotide only, news coming soon!)")
	group0.add_argument("-align", "--alignment", dest="alignment", default="", required=True, type=str, help="[MANDATORY] Input multiple sequence alignment (fasta format)")
	group0.add_argument("-o", "--output", dest="out", required=True, type=str, help="[MANDATORY] Tag for output file name")
	group0.add_argument("--site-ATCG-content", dest="N_content", required=False, default=0.0, help="[OPTIONAL: Useful to remove informative sites with excess of missing data] Minimum proportion \
						(between 0 and 1) of ATCG per informative site of the alignment (e.g. '--site-ATCG-content 1.0' will only keep positions without N's or any non-ATCG code, i.e. a core \
						SNP alignment) NOTE: This argument works on alignment positions (i.e. columns)! [default: 0.0]")
	group0.add_argument("--sample-ATCG-content", dest="ATCG_content", required=False, default=0.0, help="[OPTIONAL: Useful to remove samples with excess of missing data] Minimum proportion (\
						between 0 and 1) of ATCG in informative sites of the alignment per sample (e.g. '--ATCG-content 1.0' will only keep samples without N's or any non-ATCG code in \
						informative sites) NOTE: This argument works on samples (i.e. rows) and will be applied after '--site-ATCG-content'! [default: 0.0]")
	group0.add_argument("--remove-reference", dest="remove_ref", required=False, action="store_true", help="Set only if you want to remove the reference sequence from the alignment (reference \
						name must be provided with the argument '--reference').")
	group0.add_argument("-r", "--reference", dest="reference", required=False, type=str, default = "none", help="[OPTIONAL] Name of reference sequence. Required if '--remove-reference' and/or \
						'--use-reference-coords' specified.")
	group0.add_argument("--keep-gaps", dest="keep_gaps", required=False, action="store_true", help="Set only if you want that informative sites are kept even if they have a gap.")	
	group0.add_argument("--keep-all-gaps", dest="keep_all_gaps", required=False, action="store_true", help="Set only if you want that all sites with gaps are considered as informative.")
	group0.add_argument("--missing-code", dest="missing_code", required=False, type=str, default = "N", help="[OPTIONAL] Code representing missing data. If different from 'N', try to avoid a \
                        IUPAC character (even in lower-case) as this may influence affect the alignment cleaning. [default: N]")	
	group0.add_argument("--use-alignment-coords", dest="use_align", required=False, action="store_true", help="Set only if you want that column names in the final matrix represent the initial \
						alignment coordinates. Note: Depending on the alignment size, this argument can make alignment processing very slow!")	
	group0.add_argument("--use-reference-coords", dest="use_ref", required=False, action="store_true", help="Set only if you want that column names in the final matrix represent the reference \
						coordinates (reference name must be provided with the argument '--reference') Note: Depending on the alignment size, this argument can make alignment processing very \
						slow!")	
	group0.add_argument("-m", "--metadata", dest="metadata", required=False, default="", type=str, help="[OPTIONAL] Metadata file in .tsv format to apply sample subset.")
	group0.add_argument("-f", "--filter", dest="filter_column", required=False, default="", help="[OPTIONAL] Filter for metadata columns to select the samples of the alignment that must \
						be included in the matrix. This must be specified within quotation marks in the following format 'column< >operation< >condition' (e.g. 'country == Portugal'). When \
						more than one condition is specified for a given column, they must be separated with commas (e.g 'country == Portugal,Spain,France'). When filters include more than one \
						column, they must be separated with semicolon (e.g. 'country == Portugal,Spain,France;date > 2018-01-01;date < 2022-01-01'). White spaces are important in this argument, \
						so, do not leave spaces before and after commas/semicolons.")
	group0.add_argument("--get-positions-interest", dest="pos_int", required=False, action="store_true", help="Set only if you want a file with the final positions of interest.")
	group0.add_argument("--get-position-correspondence", dest="pos_corr", required=False, default="none", help="[OPTIONAL] Request a .tsv with position correspondence between any sequences \
						of your alignment. These should be indicated separated by a comma (e.g. seqA,seqB). To get the position coordinates of all sequences just write 'all'.")
	group0.add_argument("--position-list", dest="pos_list", required=False, default="none", help="[OPTIONAL] .tsv file with the positions of interest to be reported when \
						'--get-position-correspondence' is requested. Each column should correspond to the positions of a sequence and the sequence name should be indicated in the header. If this \
						file is not provided, all positions of the alignment will be reported.")
	group0.add_argument("--ONLY-POS-CORRESPONDENCE", dest="only_pos_corr", required=False, action="store_true", help="Set only if you JUST WANT the position correspondence and nothing else.")
	
	args = parser.parse_args()


	# starting logs	----------
	
	log_name = args.out + ".log"
	log = open(log_name, "a+")
	
	print("\n-------------------- alignment_processing.py --------------------\n")
	print("\n-------------------- alignment_processing.py --------------------\n", file = log)
	print("version", version, "last updated on", last_updated, "\n")
	print("version", version, "last updated on", last_updated, "\n", file = log)
	print(" ".join(sys.argv), "\n")
	print(" ".join(sys.argv), "\n", file = log)
	
	start = datetime.datetime.now()
	print("start:", start)
	print("start:", start, file = log)

	if args.remove_ref and args.reference == "none":
		print("'--remove-reference' specified, but I do not know what is the reference id :-(")
		print("'--remove-reference' specified, but I do not know what is the reference id :-(", file = log)
		sys.exit(1)
	elif args.use_ref and args.reference == "none":
		print("'--use-reference-coords' specified, but I do not know what is the reference id :-(")
		print("'--use-reference-coords' specified, but I do not know what is the reference id :-(", file = log)
		sys.exit(1)
		
		
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
		print("\tDetermining reference coordinates and sequence correspondence...")
		print("\tDetermining reference coordinates and sequence correspondence...", file = log)
		coords, corr_df = get_ref_coords(alignment, args.reference, args.pos_corr, args.out, args.pos_list)
		if args.only_pos_corr:
			print("You have set '--ONLY-POS-CORRESPONDENCE'. So, alignment_processing.py will exit here!")
			print("You have set '--ONLY-POS-CORRESPONDENCE'. So, alignment_processing.py will exit here!", file = log)
			sys.exit(0)
	elif args.pos_corr != "none":
		print("\tDetermining sequence correspondence...")
		print("\tDetermining sequence correspondence...", file = log)
		coords, corr_df = get_ref_coords(alignment, args.reference, args.pos_corr, args.out, args.pos_list)
		if args.only_pos_corr:
			print("You have set '--ONLY-POS-CORRESPONDENCE'. So, alignment_processing.py will exit here!")
			print("You have set '--ONLY-POS-CORRESPONDENCE'. So, alignment_processing.py will exit here!", file = log)
			sys.exit(1)

				
	# remove samples according to metadata

	if args.metadata != "" and args.filter_column != "":
		print("Filtering the alignment...")
		print("Filtering the alignment...", file = log)
		alignment = filter_align(alignment, args.metadata, args.filter_column, log)
		if len(alignment) <= 1:
			print("\nCannot proceed because " + str(len(alignment)) + " samples were kept in the matrix!")
			print("\nCannot proceed because " + str(len(alignment)) + " samples were kept in the matrix!", file = log)
			sys.exit(1)
		SeqIO.write(alignment, args.out + "_tmp.fasta", "fasta")
		alignment = AlignIO.read(args.out + "_tmp.fasta", "fasta")
		os.system("rm " + args.out + "_tmp.fasta")
	elif args.metadata != "" and args.filter_column == "":
		print("Metadata file was provided but no filter was found... I am confused :-(")
		print("Metadata file was provided but no filter was found... I am confused :-(", file = log)
		sys.exit(1)

	elif args.metadata == "" and args.filter_column != "":
		print("Metadata file was not provided but a filter was found... I am confused :-(")
		print("Metadata file was not provided but a filter was found... I am confused :-(", file = log)
		sys.exit(1)
	
	# run snp-sites to make the alignment shorter (only if '--use-reference-coords' and '--use-alignment-coords' are not set

	if float(args.ATCG_content) == 0.0 and float(args.N_content) == 1.0: # do not care about samples N content and keep only sites with ATCG
		if args.remove_ref:
			print("Removing reference sequence: " + str(args.reference))
			print("Removing reference sequence: " + str(args.reference), file = log)
			alignment = rm_ref(alignment, args.reference)
		print("Trimming the alignment with SNP-sites:")
		print("Trimming the alignment with SNP-sites:", file = log)
		print("\tsnp-sites -c " + args.out + "_tmp.fasta > " + args.out + "_tmp_flt.fasta")
		print("\tsnp-sites -c " + args.out + "_tmp.fasta > " + args.out + "_tmp_flt.fasta", file = log)
		SeqIO.write(alignment, args.out + "_tmp.fasta", "fasta")
		returned_value = os.system("snp-sites -c " + args.out + "_tmp.fasta > " + args.out + "_tmp_flt.fasta")
		if str(returned_value) != "0":
			print("\nSomething went wrong while running snp-sites :-( please double check your input files and ReporTree specifications!")
			print("\nSomething went wrong while running snp-sites :-( please double check your input files and ReporTree specifications!", file = log)
			sys.exit(1)
		returned_value = os.system("snp-sites -v -c " + args.out + "_tmp.fasta > " + args.out + "_tmp_flt.vcf")
		if str(returned_value) != "0":
			print("\nSomething went wrong while running snp-sites :-( please double check your input files and ReporTree specifications!")
			print("\nSomething went wrong while running snp-sites :-( please double check your input files and ReporTree specifications!", file = log)
			sys.exit(1)
		os.system("rm " + args.out + "_tmp.fasta")
		alignment = AlignIO.read(args.out + "_tmp_flt.fasta", "fasta")
		os.system("grep -v '##' " + args.out + "_tmp_flt.vcf > " + args.out + "_tmp.vcf")
		os.system("rm " + args.out + "_tmp_flt.fasta " + args.out + "_tmp_flt.vcf")
		print("\tAlignment length after SNP-sites: " + str(len(alignment[0].seq)))
		print("\tAlignment length after SNP-sites: " + str(len(alignment[0].seq)), file = log)

		# alignment matrix
		
		print("Getting the alignment matrix...")
		print("Getting the alignment matrix...", file = log)
		if args.use_ref:
			mx = core2mx(alignment, args.out + "_tmp.vcf", coords, log)
			os.system("rm " + args.out + "_tmp.vcf")
		else:
			mx = core2mx(alignment, args.out + "_tmp.vcf", "", log)
			os.system("rm " + args.out + "_tmp.vcf")

	else: # run snp-sites with option -c is not viable
		if args.remove_ref:
			print("Removing reference sequence: " + str(args.reference))
			print("Removing reference sequence: " + str(args.reference), file = log)
			alignment = rm_ref(alignment, args.reference)
		if not args.keep_all_gaps:
			print("Trimming the alignment with SNP-sites:")
			print("Trimming the alignment with SNP-sites:", file = log)
			print("\tsnp-sites " + args.out + "_tmp.fasta > " + args.out + "_tmp_flt.fasta")
			print("\tsnp-sites " + args.out + "_tmp.fasta > " + args.out + "_tmp_flt.fasta", file = log)
			SeqIO.write(alignment, args.out + "_tmp.fasta", "fasta")
			returned_value = os.system("snp-sites " + args.out + "_tmp.fasta > " + args.out + "_tmp_flt.fasta")
			if str(returned_value) != "0":
				print("\nSomething went wrong while running snp-sites :-( please double check your input files and ReporTree specifications!")
				print("\nSomething went wrong while running snp-sites :-( please double check your input files and ReporTree specifications!", file = log)
				sys.exit(1)
			returned_value = os.system("snp-sites -v " + args.out + "_tmp.fasta > " + args.out + "_tmp_flt.vcf")
			if str(returned_value) != "0":
				print("\nSomething went wrong while running snp-sites :-( please double check your input files and ReporTree specifications!")
				print("\nSomething went wrong while running snp-sites :-( please double check your input files and ReporTree specifications!", file = log)
				sys.exit(1)
			os.system("rm " + args.out + "_tmp.fasta")
			alignment = AlignIO.read(args.out + "_tmp_flt.fasta", "fasta")
			os.system("grep -v '##' " + args.out + "_tmp_flt.vcf > " + args.out + "_tmp.vcf")
			os.system("rm " + args.out + "_tmp_flt.fasta " + args.out + "_tmp_flt.vcf")
			print("\tAlignment length after SNP-sites: " + str(len(alignment[0].seq)))
			print("\tAlignment length after SNP-sites: " + str(len(alignment[0].seq)), file = log)

			# alignment matrix
			
			print("Getting the alignment matrix...")
			print("Getting the alignment matrix...", file = log)
			if args.use_ref:
				mx = core2mx(alignment, args.out + "_tmp.vcf", coords, log)
				os.system("rm " + args.out + "_tmp.vcf")
			else:
				mx = core2mx(alignment, args.out + "_tmp.vcf", "", log)
				os.system("rm " + args.out + "_tmp.vcf")

		else:
			# alignment matrix
			print("Getting the alignment matrix...")
			print("Getting the alignment matrix...", file = log)
			if args.use_ref:
				mx = core2mx(alignment, "", coords, log)
			else:
				mx = core2mx(alignment, "", "", log)

		# additional cleaning

		print("Clean the alignment (1st round)...")
		print("Clean the alignment (1st round)...", file = log)
		mx = clean_mx(mx, args.keep_gaps, args.keep_all_gaps, args.missing_code, log)
		run_2nd_clean = False

		# assess number of ATCG's per sample

		if float(args.ATCG_content) > 0.0:
			print("Removing samples with <= " + str(float(args.ATCG_content)*100) + "% of ATCG in informative sites...")
			print("Removing samples with <= " + str(float(args.ATCG_content)*100) + "% of ATCG in informative sites...", file = log)
			initial_samples = len(mx.index)
			mx = rm_ns(mx, args.ATCG_content, args.out, log)
			final_samples = len(mx.index)
			if len(mx.index) <= 1:
				print("Cannot proceed because " + str(len(mx.index)) + " samples were kept in the alignment!")
				print("Cannot proceed because " + str(len(mx.index)) + " samples were kept in the alignment!", file = log)
				sys.exit(1)
			if final_samples < initial_samples:
				run_2nd_clean = True

		if float(args.N_content) > 0.0 or run_2nd_clean:
			print("Clean the alignment (2nd round)...")
			print("Clean the alignment (2nd round)...", file = log)
			mx = clean_position(mx, float(args.N_content), args.keep_gaps, args.keep_all_gaps, args.missing_code, log)
			if len(mx.index) <= 1:
				print("Cannot proceed because " + str(len(mx.index)) + " samples were kept in the alignment!")
				print("Cannot proceed because " + str(len(mx.index)) + " samples were kept in the alignment!", file = log)
				sys.exit(1)
			if len(mx.columns) <= 2:
				print("Cannot proceed because " + str(len(mx.columns)-1) + " sites were kept in the alignment!")
				print("Cannot proceed because " + str(len(mx.columns)-1) + " sites were kept in the alignment!", file = log)
				sys.exit(1)	
	
	# outputs

	mx.to_csv(args.out + "_align_profile.tsv", index = False, header=True, sep ="\t")
	df2fa(mx, args.out + "_align_profile.fasta")
	print("FINAL ALIGNMENT: " + str(len(mx[mx.columns[0]].values.tolist())) + " samples and " + str(len(mx.columns) - 1) + " positions!")
	print("FINAL ALIGNMENT: " + str(len(mx[mx.columns[0]].values.tolist())) + " samples and " + str(len(mx.columns) - 1) + " positions!", file = log)
	if args.pos_int:
		pos_int(mx, args.out + "_positions_of_interest.tsv")

end = datetime.datetime.now()
elapsed = end - start

print("\nalignment_processing.py is done!")
print("\nalignment_processing.py is done!", file = log)
print("\nEnd:", end)
print("\nEnd:", end, file = log)
print("Time elapsed:", elapsed)
print("Time elapsed:", elapsed, file = log)

log.close()
