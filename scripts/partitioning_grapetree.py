#!/usr/bin/env	python3

"""
Obtain genetic clusters at any partition level(s) of a minimum spanning tree derived from a cg/wgMLST allele matrix.
It requires a MODIFIED version of GrapeTree available at https://github.com/insapathogenomics/GrapeTree
Note: This script takes advantage of GrapeTree -> do not forget to cite its authors as well!!
By Veronica Mixao
@INSA
"""

import sys
import os
import argparse
import textwrap
import glob
import pandas
from datetime import date
import datetime as datetime
import string

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

partitioning_grapetree_script = os.path.realpath(__file__)
grapetree = partitioning_grapetree_script.rsplit("/", 1)[0] + "/GrapeTree/grapetree.py"
python = sys.executable

version = "1.5.0"
last_updated = "2024-07-12"

# additional functions	----------

def get_loci2use(loci,allele_profile):
	""" get the list of loci to include """
	
	mx = allele_profile
	loci2include = []
	with open(loci) as inloci:
		lines = inloci.readlines()
		for line in lines:
			l =line.split("\n")[0]
			if l not in loci2include:
				if l in mx.columns:
					loci2include.append(l)
				else:
					sys.exit("Locus " + str(l) + " is not present in the provided allele matrix. Please revise your list before we proceed!!")
			else:
				sys.exit("Locus " + str(l) + " is duplicated in the loci list. Please revise your list before we proceed!!")
	
	return loci2include


# running the pipeline	----------

def main():
	# defining parameters ----------

	parser = argparse.ArgumentParser(prog="partitioning_grapetree.py", formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent("""\
										###############################################################################             
										#                                                                             #
										#                          partitioning_grapetree.py                          #
										#                                                                             #
										###############################################################################
																	
										partitioning_grapetree.py obtains genetic clusters at any partition level(s) of 
										a minimum spanning tree derived from a SNP or cg/wgMLST allele matrix.
										
										
										This script requires a modified version of GrapeTree that is available at: 
										https://github.com/insapathogenomics/GrapeTree
										
										
										NOTE: Do not forget to also cite GrapeTree original authors.
										
										
										How to run partitioning_grapetree.py?
										
										A) Minimum-spanning tree of the provided allele matrix and partitions for 
										thresholds 5,8, and 10 to 20:
										partitioning_grapetree.py -a ALLELE_PROFILE -o OUTPUT_NAME --method MSTreeV2 
										--missing 0 --n_proc 5 -thr 5,8,10-20 
										
										B) Minimum-spanning tree for a subset of the provided allele matrix and all 
										partitions:
										partitioning_grapetree.py -a ALLELE_PROFILE -o OUTPUT_NAME --method MSTreeV2 
										--missing 0 --n_proc 5 -thr max -m METADATA -f 
										"column_metadata<>operation<>value"
										
										-------------------------------------------------------------------------------"""))

	group0 = parser.add_argument_group("Partitioning with GrapeTree", "Specifications to get and cut minimum spanning trees")
	group0.add_argument("-a", "--allele-profile", dest="allele_profile", required=True, type=str, help="[MANDATORY] Input profile matrix (can either be an allele matrix or a SNP matrix)")
	group0.add_argument("-l", "--loci", dest="loci", required=False, type=str, default = "none", help="[OPTIONAL] List of loci (e.g. cgMLST) that must be used for the clustering analysis. If \
					 	'--site-inclusion' argument > 0, this list of loci will be complemented with additional loci that fulfill the requirements specified in this argument.")
	group0.add_argument("-o", "--output", dest="out", required=True, type=str, help="[MANDATORY] Tag for output file name")
	group0.add_argument("--site-inclusion", dest="samples_called", required=False, default = 0.0, help="[OPTIONAL: Useful to remove informative sites/loci with excess of missing data] Minimum \
						proportion of samples per site/loci without missing data (e.g. '--site-inclusion 1.0' will only keep loci/positions without missing data, i.e. a core alignment; \
						'--site-inclusion 0.0' will keep all loci/positions) NOTE: This argument works on profile/alignment loci/positions (i.e. columns)! [default: 0.0].")
	group0.add_argument("--loci-called", dest="loci_called", required=False, default = 0.0, help="[OPTIONAL] Minimum percentage of loci/positions called for allele/SNP matrices (e.g. \
						'--loci-called 0.95' will only keep in the profile matrix samples with > 95%% of alleles/positions, i.e. <= 5%% missing data). Applied after '--site-inclusion' argument! \
						[default: 0.0]")
	group0.add_argument("--method", dest="grapetree_method", default="MSTreeV2", help="\"MSTreeV2\" [DEFAULT]\n Alternative:\"MSTree\" (goeBURST)\n(Note: When using profile without missing data, \
						we suggest that you use MSTree)")
	group0.add_argument("--missing", dest="handler", default=0, type=int, help="ONLY FOR MSTree. \n0: [DEFAULT] ignore missing data in pairwise comparison. \n1: remove column \
						with missing data. \n2: treat missing data as an allele. \n3: use absolute number of allelic differences.")
	group0.add_argument("--missing-code", dest="missing_code", required=False, type=str, default = "0", help="[OPTIONAL] Code representing missing data. If different from '0' or 'N', please try \
						to avoid a IUPAC character (even in lower-case). [default: 0]")
	group0.add_argument("--wgMLST", dest="wgmlst", default=False, action="store_true", help="[EXPERIMENTAL] a better support of wgMLST schemes (check GrapeTree github for details).")
	group0.add_argument("--n_proc",  dest="number_of_processes", type=int, default=5, help="Number of CPU processes in parallel use. [5]")
	group0.add_argument("-thr", "--threshold", dest="threshold", default = "max", help="[OPTIONAL] Partition thresholds for clustering definition. Different thresholds can be comma-separated \
						(e.g. 5,8,16). Ranges can be specified with an hyphen separating minimum and maximum (e.g. 5,8,10-20). If this option is not set, the script will perform clustering for \
						all the values in the range 0 to max. Note: Threshold values are inclusive, i.e. '-thr 7' will consider samples with <= 7 differences as belonging to the same cluster!")
	group0.add_argument("-pct_thr", "--pct_threshold", dest="pct_threshold", default = "none", help="[[OPTIONAL] Similar to 'thr' but values are indicated as the proportion of differences to the \
						final allelic schema size or number of informative positions, e.g. '-pct_thr 0.005' corresponds to a threshold of 5 allelic/SNP differences in a matrix with 1000 \
						loci/sites under analysis). Different values can be comma-separated (e.g. 0.005,0.01,0.1). Ranges CANNOT be specified. This option is particularly useful for dynamic \
						wgMLST analysis for which the size of the schema under analysis is contigent on dataset diversity. Note: This argument can be specified even if you used the '-thr' \
						argument.")
	group0.add_argument("-m", "--metadata", dest="metadata", required=False, default="", type=str, help="[OPTIONAL] Metadata file in .tsv format to select the samples to reconstruct the minimum \
						spanning tree according to the '--filter' argument")
	group0.add_argument("-f", "--filter", dest="filter_column", required=False, default="", help="[OPTIONAL] Filter for metadata columns to select the samples of the allele matrix that must \
						be used for tree reconstruction. This must be specified within quotation marks in the following format 'column< >operation< >condition' (e.g. 'country == Portugal'). When \
						more than one condition is specified for a given column, they must be separated with commas (e.g 'country == Portugal,Spain,France'). When filters include more than one \
						column, they must be separated with semicolon (e.g. 'country == Portugal,Spain,France;date > 2018-01-01;date < 2022-01-01'). White spaces are important in this argument, \
						so, do not leave spaces before and after commas/semicolons.")
	group0.add_argument("-d", "--dist", dest="dist", required=False, default=1.0, type=float, help="Distance unit by which partition thresholds will be multiplied (example: if -d 10 and \
						-thr 5,8,10-30, the tree will be cut at 50,80,100,110,120,...,300). Currently, the default is 1, which is equivalent to 1 allele/SNP distance. [1.0]")
	group0.add_argument("--matrix-4-grapetree", dest="matrix4grapetree", required=False, action="store_true", help="Output an additional allele profile matrix with the header ready for GrapeTree \
						visualization. Set only if you WANT the file!")	
							
	args = parser.parse_args()


	# starting logs	----------

	log_name = args.out + ".log"
	log = open(log_name, "a+")
		
	print("\n-------------------- partitioning_grapetree.py --------------------\n")
	print("\n-------------------- partitioning_grapetree.py --------------------\n", file = log)
	print("version", version, "last updated on", last_updated, "\n")
	print("version", version, "last updated on", last_updated, "\n", file = log)
	print(" ".join(sys.argv), "\n")
	print(" ".join(sys.argv), "\n", file = log)

	start = datetime.datetime.now()
	print("start:", start)
	print("start:", start, file = log)

	missing_code = str(args.missing_code)
	if missing_code != "0":
		missing_need = True
	else:
		missing_need = False

	# filtering allele matrix	----------

	if args.metadata != "" and args.filter_column != "" and ".fasta" not in args.allele_profile and ".fa" not in args.allele_profile and ".fas" not in args.allele_profile:
		print("Filtering the allele matrix...")
		print("Filtering the allele matrix...", file = log)
		
		filters = args.filter_column
		mx = pandas.read_table(args.metadata, dtype = str)
		columns_names = [col.replace(" ", "_") for col in mx.columns]
		mx.replace(" ", "_", regex=True, inplace=True)
		mx.columns = columns_names
		allele_mx = pandas.read_table(args.allele_profile, dtype = str)
		allele_mx = allele_mx.replace("INF-","", regex=True) #implemented because of chewie-ns profiles
		allele_mx = allele_mx.replace("\*","", regex=True) #implemented because of chewie-ns profiles
		if missing_need:
			allele_mx_id = allele_mx[allele_mx.columns[0]]
			if missing_code == "empty":
				allele_mx.fillna("0", inplace = True)
			else:
				allele_mx = allele_mx.replace({missing_code: "0"})
			allele_mx[allele_mx.columns[0]] = allele_mx_id
			allele_mx.to_csv(args.out + "_missing0.tsv", index = False, header=True, sep ="\t")
			allele_filename = args.out + "_missing0.tsv"
			missing_need = False
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
						print("\tCannot apply the 'iso_week' filter because the column 'date' was not found in the metadata!")
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
		
		initial_samples = len(allele_mx[allele_mx.columns[0]].values.tolist())
		allele_mx_filtered = allele_mx[allele_mx[allele_mx.columns[0]].isin(samples)]
		final_samples = len(allele_mx_filtered[allele_mx_filtered.columns[0]].values.tolist())
		print("\tFrom the " + str(initial_samples) + " samples, " + str(final_samples) + " were kept in the matrix...")
		print("\tFrom the " + str(initial_samples) + " samples, " + str(final_samples) + " were kept in the matrix...", file = log)
		if final_samples <= 1:
			print("\nCannot proceed because " + str(final_samples) + " samples were kept in the matrix!")
			print("\nCannot proceed because " + str(final_samples) + " samples were kept in the matrix!", file = log)
			sys.exit(1)
		allele_mx_filtered.to_csv(args.out + "_subset_matrix.tsv", index = False, header=True, sep ="\t")
		allele_filename = args.out + "_subset_matrix.tsv"
		total_size = len(allele_mx_filtered.columns) - 1
		
	elif args.metadata != "" and args.filter_column == "":
		print("Metadata file was provided but no filter was found... I am confused :-(")
		print("Metadata file was provided but no filter was found... I am confused :-(", file = log)
		sys.exit()

	elif args.metadata == "" and args.filter_column != "":
		print("Metadata file was not provided but a filter was found... I am confused :-(")
		print("Metadata file was not provided but a filter was found... I am confused :-(", file = log)
		sys.exit()

	else:
		sample_column = "sequence"
		allele_filename = args.allele_profile
		allele_mx = pandas.read_table(allele_filename, dtype = str)
		allele_mx = allele_mx.replace("INF-","", regex=True) #implemented because of chewie-ns profiles
		allele_mx = allele_mx.replace("\*","", regex=True) #implemented because of chewie-ns profiles
		allele_mx.to_csv(args.out + "_temporary_clean_codes.tsv", index = False, header=True, sep ="\t")
		allele_filename = args.out + "_temporary_clean_codes.tsv"
		if missing_need:
			allele_mx_id = allele_mx[allele_mx.columns[0]]
			if missing_code == "empty":
				allele_mx.fillna("0", inplace = True)
			else:
				allele_mx = allele_mx.replace({str(missing_code): "0"})
			allele_mx[allele_mx.columns[0]] = allele_mx_id
			allele_mx.to_csv(args.out + "_missing0.tsv", index = False, header=True, sep ="\t")
			allele_filename = args.out + "_missing0.tsv"
			missing_need = False
		total_size = len(allele_mx.columns) - 1


	# cleaning allele matrix (columns)
	mx_allele = pandas.read_table(allele_filename, dtype = str)
	pos_t0 = len(mx_allele.columns[1:])
	if args.loci != "none":
		loci2include = get_loci2use(args.loci,mx_allele)
		if float(args.samples_called) == 0.0:
			print("Keeping the sites/loci present in the loci list.")
			print("Keeping the sites/loci present in the loci list.", file = log)
			for col in mx_allele.columns[1:]:
				if col not in loci2include:
					mx_allele = mx_allele.drop(columns=col)
		else:
			print("Keeping the sites/loci present in the loci list and those with information in >= " + str(float(args.samples_called) * 100) + "%% of the samples...")
			print("Keeping the sites/loci present in the loci list and those with information in >= " + str(float(args.samples_called) * 100) + "%% of the samples...", file = log)
			for col in mx_allele.columns[1:]:
				if col not in loci2include:
					values = mx_allele[col].values.tolist()
					if (len(values)-values.count("0"))/len(values) < float(args.samples_called):
						mx_allele = mx_allele.drop(columns=col)
					elif (len(values)-values.count(0))/len(values) < float(args.samples_called):
						mx_allele = mx_allele.drop(columns=col)
		pos_t1 = len(mx_allele.columns[1:])
		print("\tFrom the " + str(pos_t0) + " loci/positions, " + str(pos_t1) + " were kept in the matrix.")
		print("\tFrom the " + str(pos_t0) + " loci/positions, " + str(pos_t1) + " were kept in the matrix.", file = log)
	else:
		pos_t1 = len(mx_allele.columns[1:])
		if float(args.samples_called) > 0.0:
			print("Keeping the sites/loci with information in >= " + str(float(args.samples_called) * 100) + "%% of the samples...")
			print("Keeping the sites/loci with information in >= " + str(float(args.samples_called) * 100) + "%% of the samples...", file = log)
			for col in mx_allele.columns[1:]:
				values = mx_allele[col].values.tolist()
				if (len(values)-values.count("0"))/len(values) < float(args.samples_called):
					mx_allele = mx_allele.drop(columns=col)
				elif (len(values)-values.count(0))/len(values) < float(args.samples_called):
					mx_allele = mx_allele.drop(columns=col)	
			pos_t1 = len(mx_allele.columns[1:])
			print("\tFrom the " + str(pos_t0) + " loci/positions, " + str(pos_t1) + " were kept in the matrix.")
			print("\tFrom the " + str(pos_t0) + " loci/positions, " + str(pos_t1) + " were kept in the matrix.", file = log)
		
	if pos_t1 <= 1:
		print("\nCannot proceed because " + str(pos_t1) + " positions/alleles were kept in the matrix!")
		print("\nCannot proceed because " + str(pos_t1) + " positions/alleles were kept in the matrix!", file = log)
		sys.exit()
	with open(args.out + "_loci_used.txt", "w+") as loci_out:
		for locus in mx_allele.columns[1:]:
			print(locus, file = loci_out)
	mx_allele.to_csv(args.out + "_clean_missing_matrix.tsv", index = False, header=True, sep ="\t")
	allele_filename = args.out + "_clean_missing_matrix.tsv"
	total_size = len(mx_allele.columns) - 1
			

	# cleaning allele matrix (rows)	----------

	if float(args.loci_called) > 0.0:
		print("Cleaning the profile matrix using a threshold of >" + str(args.loci_called) + " alleles/positions called...")
		print("Cleaning the profile matrix using a threshold of >" + str(args.loci_called) + " alleles/positions called...", file = log)
		
		allele_mx = pandas.read_table(allele_filename, dtype = str)
		report_allele_mx = {}
		
		len_schema = len(allele_mx.columns) - 1
		
		report_allele_mx["samples"] = allele_mx[allele_mx.columns[0]]
		report_allele_mx["missing"] = allele_mx.isin(["0"]).sum(axis=1)
		report_allele_mx["called"] = len_schema - allele_mx.isin(["0"]).sum(axis=1)
		report_allele_mx["pct_called"] = (len_schema - allele_mx.isin(["0"]).sum(axis=1)) / len_schema

		report_allele_df = pandas.DataFrame(data = report_allele_mx)
		if float(args.loci_called) != 1.0:
			flt_report = report_allele_df[report_allele_df["pct_called"] > float(args.loci_called)]
		else:
			flt_report = report_allele_df[report_allele_df["pct_called"] == float(args.loci_called)]
		pass_samples = flt_report["samples"].values.tolist()
		
		print("\tFrom the " + str(len(allele_mx[allele_mx.columns[0]].values.tolist())) + " samples, " + str(len(pass_samples)) + " were kept in the profile matrix.")
		print("\tFrom the " + str(len(allele_mx[allele_mx.columns[0]].values.tolist())) + " samples, " + str(len(pass_samples)) + " were kept in the profile matrix.", file = log)
		
		if len(pass_samples) <= 1:
			print("\nCannot proceed because " + str(len(pass_samples)) + " samples were kept in the matrix!")
			print("\nCannot proceed because " + str(len(pass_samples)) + " samples were kept in the matrix!", file = log)
			sys.exit()
			
		allele_mx = allele_mx[allele_mx[allele_mx.columns[0]].isin(pass_samples)]
		allele_mx.to_csv(args.out + "_flt_samples_matrix.tsv", index = False, header=True, sep ="\t")
		allele_filename = args.out + "_flt_samples_matrix.tsv"
		report_allele_df.to_csv(args.out + "_loci_report.tsv", index = False, header=True, sep ="\t")
		total_size = len(allele_mx.columns) - 1
		
	# preparing allele matrix for grapetree	----------

	if args.matrix4grapetree:
		mx_allele = pandas.read_table(allele_filename, dtype = str)
		first_col = str(mx_allele.columns[0])
		if first_col[0] != "#":
			mx_allele = mx_allele.rename(columns={first_col: "#" + first_col})
			mx_allele.to_csv(args.out + "_alleles_4_grapetree.tsv", index = False, header=True, sep ="\t")

	
	# getting distance matrix - hamming	----------

	print("Getting the pairwise distance matrix with cgmlst-dists (if your profile matrix is too big, this will be done in chunks of 2000 alleles/positions)...")
	print("Getting the pairwise distance matrix with cgmlst-dists (if your profile matrix is too big, this will be done in chunks of 2000 alleles/positions)...", file = log)
							
	# convert ATCG to integers
	allele_mx = pandas.read_table(allele_filename, dtype = str)
	initial_ids = allele_mx[allele_mx.columns[0]]

	string_values_lower = list(string.ascii_lowercase)
	string_values_upper = list(string.ascii_uppercase)
	string_values = string_values_lower + string_values_upper

	str2int = {}
	count_string = 1
	for info in string_values:
		str2int[info] = count_string
		count_string += 1
	if missing_need:
		str2int[missing_code] = "0"

	alleles = allele_mx.replace(str2int)
	alleles[alleles.columns[0]] = initial_ids

	# divide a big dataframe into chunks
	df_size = len(alleles.columns) - 1 
	start_chunk = 0
	chunk_size = 2000
	end_chunk = start_chunk + chunk_size
	df_counter = 0

	while end_chunk < df_size + chunk_size:
		main_df = alleles.set_index(alleles.columns[0], drop = True)
		tmp_df = main_df.iloc[: , start_chunk:end_chunk]
		df_counter += 1
		tmp_df.to_csv(args.out + "_temporary_profile.tsv", index = True, header = True, sep ="\t")

		# run cgmlst-dists
		returned_value = os.system("cgmlst-dists " + args.out + "_temporary_profile.tsv > " + args.out + "_tmp_dist_hamming.tsv")
		if str(returned_value) != "0":
			print("\nSomething went wrong while running cgmlst-dists to get hamming distances :-( please double check your input files and ReporTree specifications!")
			print("\nSomething went wrong while running cgmlst-dists to get hamming distances :-( please double check your input files and ReporTree specifications!", file = log)
			sys.exit(1)
		os.system("rm " + args.out + "_temporary_profile.tsv")
		sub_dist_df = pandas.read_table(args.out + "_tmp_dist_hamming.tsv")
		os.system("rm " + args.out + "_tmp_dist_hamming.tsv")
		sub_dist_df.set_index(sub_dist_df.columns[0], inplace = True, drop = True)
		if df_counter == 1:
			dist_df = sub_dist_df
		else:
			dist_df = dist_df.add(sub_dist_df, axis = 0)
		start_chunk += chunk_size
		end_chunk += chunk_size

	dist_df.index.names = ["dists"]	
	dist_df.to_csv(args.out + "_dist_hamming.tsv", sep = "\t", index = True, header = True)					
	

	# running grapetree	----------

	print("Running GrapeTree...")
	print("Running GrapeTree...", file = log)

	extra_commands = ""
	if args.wgmlst == True:
		extra_commands += " --wgMLST" 

	print(python + " " + grapetree + " -p " + allele_filename + " -m " + args.grapetree_method + " -o " + args.out + " --missing " + str(args.handler) + " --n_proc " + str(args.number_of_processes) + extra_commands)
	print(python + " " + grapetree + " -p " + allele_filename + " -m " + args.grapetree_method + " -o " + args.out + " --missing " + str(args.handler) + " --n_proc " + str(args.number_of_processes) + extra_commands, file = log)
	returned_value = os.system(python + " " + grapetree + " -p " + allele_filename + " -m " + args.grapetree_method + " -o " + args.out + " --missing " + str(args.handler) + " --n_proc " + str(args.number_of_processes) + extra_commands)
	if str(returned_value) != "0":
		print("\nSomething went wrong while running GrapeTree :-( Please check your input matrix and your filtering options... ")
		print("\nSomething went wrong while running GrapeTree :-( Please check your input matrix and your filtering options... ", file = log)
		sys.exit(1)


	# defining the threshold range ----------

	print("\nProcessing clustering threshold...")
	print("\nProcessing clustering threshold...", file = log)

	distances = []
	range_def = []

	with open(args.out + "_dist.tsv", "r") as groups:
		group = groups.readlines()
		for line in group:
			lin = line.split("\n")
			s1,s2,d = lin[0].split("\t")
			distances.append(int(d))
		min_distance = 0
		max_distance = max(distances) + 1
			
	distance_magnitude = args.dist


	## look at specific thresholds

	if args.threshold != "max":
		if "," in args.threshold:
			values = args.threshold.split(",")
			for v in values:
				if "-" in v:
					min_range = int(v.split("-")[0])
					max_range = int(v.split("-")[1]) + 1
					r = min_range,max_range
					range_def.append(r)
				else:
					min_range = int(v)
					max_range = int(v) + 1
					r = min_range,max_range
					range_def.append(r)
		else:
			v = args.threshold
			if "-" in v:
				min_range = int(v.split("-")[0])
				max_range = int(v.split("-")[1]) + 1
				r = min_range,max_range
				range_def.append(r)
			else:
				v = args.threshold
				min_range = int(v)
				max_range = int(v) + 1
				r = min_range,max_range
				range_def.append(r)

	else:
		r = min_distance,max_distance
		range_def.append(r)


	## look at pct thresholds

	range_def_pct = []
	pct_correspondence = {}
	if args.pct_threshold != "none":
		print("\n\tCorrespondence between percentage and number of differences:")
		print("\n\tCorrespondence between percentage and number of differences:", file = log)
		print("\n\t#PERCENTAGE\tDIFFERENCES")
		print("\n\t#PERCENTAGE\tDIFFERENCES", file = log)
		pcts = args.pct_threshold.split(",")
		for pct in pcts:
			pct_threshold = int(float(pct)*float(total_size))
			print("\t" + str(float(pct)*100) + "\t" + str(pct_threshold))
			print("\t" + str(float(pct)*100) + "\t" + str(pct_threshold), file = log)
			min_range = pct_threshold
			max_range = pct_threshold + 1
			r = min_range,max_range		
			if pct_threshold not in pct_correspondence.keys():
				pct_correspondence[pct_threshold] = []
				range_def_pct.append(r)
			pct_correspondence[pct_threshold].append(str(pct))

	all_thresholds = {"general" : range_def, "pct" : range_def_pct}
			
		
	# processing redundant info ----------

	print("Getting redundant sample information...")
	print("Getting redundant sample information...", file = log)
	red = {} #dictionary with redundant names

	with open(args.out + "_redundantSamples.txt", "r") as redundant_file:
		redundant = redundant_file.readlines()
		
		for line in redundant:
			l = line.split("\n")[0]
			representative = l.split("\t")[0]
			others = l.split("\t")[1]
			additional = []
			if "," in others:
				for sample in others.split(","):
					additional.append(sample)			
				red[representative] = additional


	# defining clusters ----------

	print("Defining clusters...")
	print("Defining clusters...", file = log)
	info = {} #dictionary with all the partitions -> info[partition][cluster] = set(sample1, sample2,...)
	order_partitions = [] #list of partitions... useful to write the dataframe with pandas
	order_partitions.append(sample_column)

	for thr_type in all_thresholds.keys():
		range_list = all_thresholds[thr_type]
		if len(range_list) > 0:
			for min_r,max_r in range_list:
				if thr_type == "general":
					print("\tCalculating clustering in range",min_r,max_r,"with a distance of",distance_magnitude)
					print("\tCalculating clustering in range",min_r,max_r,"with a distance of",distance_magnitude, file = log)
				elif thr_type == "pct":
					print("\tCalculating clustering for partition " + str(min_r) + ", which corresponds to the pct threshold of: " + ", ".join(pct_correspondence[int(min_r)]))
					print("\tCalculating clustering for partition " + str(min_r) + ", which corresponds to the pct threshold of: " + ", ".join(pct_correspondence[int(min_r)]), file = log)
				for partition_value in range(min_r,max_r):
					if thr_type == "general":
						partition = partition_value * distance_magnitude
					elif thr_type == "pct":
						partition = partition_value
					if partition >= max_distance:
						print("\t\tThe requested partition, " + str(int(partition)) +  ", is higher than the maximum distance. Clustering will not be computed.")
						print("\t\tThe requested partition, " + str(int(partition)) +  ", is higher than the maximum distance. Clustering will not be computed.", file = log)
						break
					else:
						if thr_type == "general":
							order_partitions.append("MST-" + str(partition_value) + "x" + str(distance_magnitude))
						elif thr_type == "pct":
							for percentage_thresholds in pct_correspondence[int(min_r)]:
								order_partitions.append("MST-" + str(percentage_thresholds))
						clusters = {} #dictionary with the different clusters -> clusters[cluster] = set(sample1, sample2,...)
						i = 0 #cluster name definer
						checked_samples = set() #samples that have been checked
						
						with open(args.out + "_dist.tsv", "r") as group_file:
							group = group_file.readlines()
							for line in group:
								lin = line.split("\n")
								s1,s2,d = lin[0].split("\t") #sample1,sample2,distance

								if int(d) <= partition: #they are a cluster
									if s1 not in checked_samples and s2 not in checked_samples: #s1,s2 do not have an assigned cluster and should be added to clusters dictionary with same key
										i += 1
										clusters[i] = set()
										clusters[i].add(s1)
										clusters[i].add(s2)
										checked_samples.add(s1)
										checked_samples.add(s2)
									else:
										if s1 in checked_samples and s2 not in checked_samples: #s1 has been assigned to a cluster before
											for cluster in clusters.keys():
												if s1 in clusters[cluster]: #finding the cluster where s1 belongs
													clusters[cluster].add(s2) #add s2
													checked_samples.add(s2)
										else:
											if s2 in checked_samples and s1 not in checked_samples: #s2 has been assigned to a cluster before
												for cluster in clusters.keys():
													if s2 in clusters[cluster]: #finding the cluster where s2 belongs
														clusters[cluster].add(s1) #add s1
														checked_samples.add(s1)
											else: 
												if s1 in checked_samples and s2 in checked_samples: #both samples have been assigned to a cluster before
													for cluster in clusters.keys():
														if s1 in clusters[cluster]:
															cluster_s1 = cluster
														if s2 in clusters[cluster]:
															cluster_s2 = cluster
													if cluster_s1 != cluster_s2: #was it the same cluster? 
														i += 1
														clusters[i] = clusters[cluster_s1] | clusters[cluster_s2] #union of the two clusters
														del clusters[cluster_s1] #delete previous clusters
														del clusters[cluster_s2] #delete previous clusters
								else: #they are not a cluster
									if s1 not in checked_samples and s2 not in checked_samples: #they have not a cluster assigned yet
										i += 1
										clusters[i] = set()
										clusters[i].add(s1) #create cluster for s1
										i += 1
										clusters[i] = set()
										clusters[i].add(s2) #create cluster for s2
										checked_samples.add(s1)
										checked_samples.add(s2)
									else:
										if s1 in checked_samples and s2 not in checked_samples: #s1 has been assigned to a cluster before
											i += 1
											clusters[i] = set()
											clusters[i].add(s2) #create a cluster just for s2
											checked_samples.add(s2)
										else:
											if s2 in checked_samples and s1 not in checked_samples: #s2 has been assigned to a cluster before
												i += 1
												clusters[i] = set()
												clusters[i].add(s1) #create a cluster just for s2
												checked_samples.add(s1)
											else:
												if s1 in checked_samples and s2 in checked_samples: #both samples have been assigned to a cluster before
													for cluster in clusters.keys():
														if s1 in clusters[cluster]:
															cluster_s1 = cluster
														if s2 in clusters[cluster]:
															cluster_s2 = cluster
					if thr_type == "general":
						info["MST-" + str(partition_value) + "x" + str(distance_magnitude)] = clusters
					elif thr_type == "pct":
						for percentage_thresholds in pct_correspondence[int(min_r)]:
							info["MST-" + str(percentage_thresholds)] = clusters


	info_sample = {} #dictionary with all the samples -> info_sample[sample][partition] = cluster

	with open(args.out + "_clusterComposition.tsv", "w+") as clusterComposition:
		print("Creating cluster composition file...")
		print("Creating cluster composition file...", file = log)
		print("#partition\tcluster\tcluster_length\tsamples", file = clusterComposition)
		for partition in info.keys():
			cluster_counter = 0 #to rename cluster
			singleton_counter = 0 #to rename and distinguish singletons from clusters
			for cluster in info[partition].keys():
				cluster_composition = []
				for sample in info[partition][cluster]:
					if sample in red.keys(): #if it has redundants
						for red_sample in red[sample]:
							cluster_composition.append(red_sample)		
					else:
						cluster_composition.append(sample)

				if len(cluster_composition) == 0:
					print("\tERROR!! CLUSTER COMPOSITION == 0!!!")
				elif len(cluster_composition) == 1:
					singleton_counter += 1
					name = "singleton_" + str(singleton_counter)
				else:
					cluster_counter += 1
					name = "cluster_" + str(cluster_counter)
					
				print(str(partition) + "\t" + str(name) + "\t" + str(len(cluster_composition)) + "\t" + ",".join(cluster_composition), file = clusterComposition)

				for s in cluster_composition:
					if s not in info_sample.keys():
						info_sample[s] = {}
					info_sample[s][partition] = name
					

	# getting partitions file ----------

	print("Creating sample partitions file...")
	print("Creating sample partitions file...", file = log)

	typing = {}
	typing[sample_column] = []

	for sample in info_sample.keys():
		typing[sample_column].append(sample)
		for thr_type in all_thresholds.keys():
			range_list = all_thresholds[thr_type]
			for min_r,max_r in range_list:
				for partition_value in range(min_r,max_r):
					if partition_value < max_distance:
						if thr_type == "general":
							partition = partition_value * distance_magnitude
							if "MST-" + str(partition_value) + "x" + str(distance_magnitude) not in typing.keys():
								typing["MST-" + str(partition_value) + "x" + str(distance_magnitude)] = []
							typing["MST-" + str(partition_value) + "x" + str(distance_magnitude)].append(info_sample[sample]["MST-" + str(partition_value) + "x" + str(distance_magnitude)])
						elif thr_type == "pct":
							partition = partition_value
							for percentage_thresholds in pct_correspondence[int(partition)]:
								if "MST-" + str(percentage_thresholds) not in typing.keys():
									typing["MST-" + str(percentage_thresholds)] = []
								typing["MST-" + str(percentage_thresholds)].append(info_sample[sample]["MST-" + str(percentage_thresholds)])


	#output percentage correspondence

	if len(pct_correspondence.keys()) > 0:
		with open(args.out + "_parcentage2threshold.tsv", "w+") as out_pct:
			print("#percentage\tthreshold", file = out_pct)
			for threshold in pct_correspondence.keys():
				for percentage in pct_correspondence[threshold]:
					print(str(percentage) + "\t" + str(threshold), file = out_pct)

	if os.path.exists(args.out + "_temporary_clean_codes.tsv"):
		os.system("rm " + args.out + "_temporary_clean_codes.tsv")
		
	matrix = pandas.DataFrame(data = typing, columns = order_partitions)
	matrix.to_csv(args.out + "_partitions.tsv", index = False, header=True, sep ="\t")

	end = datetime.datetime.now()
	elapsed = end - start

	print("\npartitioning_grapetree.py is done!")
	print("\npartitioning_grapetree.py is done!", file = log)
	print("\nEnd:", end)
	print("\nEnd:", end, file = log)
	print("Time elapsed:", elapsed)
	print("Time elapsed:", elapsed, file = log)

	log.close()

if __name__ == "__main__":
    main()
