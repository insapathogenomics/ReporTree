#!/usr/bin/env	python3

"""
Obtain genetic clusters at any partition level(s) of a distance matrixes using hierarchical clustering
By Veronica Mixao
@INSA
"""

import sys
import os
import argparse
import textwrap
import pandas
from datetime import date
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, maxdists, to_tree
from scipy.spatial.distance import squareform
import string

partitioning_HC_script = os.path.realpath(__file__)

sys.setrecursionlimit(10000) # please increase this number, if you are getting the error "RecursionError: maximum recursion depth exceeded while calling a Python object" 

version = "1.7.0"
last_updated = "2024-07-12"

# functions	----------

def conv_nucl(alleles, missing_code, missing_need):
	"""convert nucl to integers"""
	
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
	alleles = alleles.replace(str2int)

	return alleles
	
	
def filter_mx(matrix, mx, filters, matrix_type, log):
	""" filter the allele or pairwise distance matrix
	input: matrix
	output: filtered pandas dataframe
	"""
    
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
	if matrix_type == "dist":
		pairwise_dist_mx = matrix.set_index(matrix.columns[0], drop = True)
		columns_interest = []
		
		for sample in pairwise_dist_mx.columns:
			if sample in samples:
				columns_interest.append(sample)

		df = pairwise_dist_mx[columns_interest].loc[columns_interest]
		df = df.reset_index(drop=False)
	
	else:
		df = matrix[matrix[matrix.columns[0]].isin(samples)]
	
	return df

			
def dist_mx(dist, log):
    """ obtain a condensed distance matrix
    input: squared pairwise distance matrix
    output: consensed distance matrix
    """
    
    print("\tGetting condensed distance matrix...")
    print("\tGetting condensed distance matrix...", file = log)

    dist = dist.set_index(dist.columns[0],drop=True)
    
    samples = dist.columns
    
    condensed_dist = squareform(dist)
    
    return dist, condensed_dist, samples
    
    			
def hcluster(dist_mx, method_choice, log):
    """ obtain linkage array 
    input: distance matrix
    output: linkage array, list samples names, max dist
    """
    
    clustering = linkage(dist_mx, method = method_choice)
    max_dist = maxdists(clustering)[-1]
    
    print("\tMaximum distance " + str(max_dist) + "...")
    print("\tMaximum distance " + str(max_dist) + "...", file = log)
    
    return clustering, max_dist
	

def get_partitions(clustering, threshold, log):    
    """ obtain clustering information for a given threshold
    input: linkage array
    output: list with cluster names
    """
    
    clusters = list(fcluster(clustering, t = threshold, criterion = "distance"))
    
    return clusters


def get_cluster_composition(outfile, partitions):
	""" get summary of cluster composition 
	input: clusters dataframe
	output: summary dataframe 
	"""
	
	with open(outfile, "w+") as output:
		print("#partition\tcluster\tcluster_length\tsamples", file = output)
		for partition in partitions.keys():
			for cluster in partitions[partition].keys():
				for cluster_length in partitions[partition][cluster].keys():
					print(partition + "\t" + cluster + "\t" + str(cluster_length) + "\t" + ",".join(partitions[partition][cluster][cluster_length]), file = output)

	return 	

def get_newick(node, parent_dist, leaf_names, newick='') -> str:
    """
    This function was retrieved from: https://stackoverflow.com/questions/28222179/save-dendrogram-to-newick-format
    We thank STACKOVERFLOW and @MrTomRod and @jfn
    
    Convert sciply.cluster.hierarchy.to_tree()-output to Newick format.
    :param node: output of sciply.cluster.hierarchy.to_tree()
    :param parent_dist: output of sciply.cluster.hierarchy.to_tree().dist
    :param leaf_names: list of leaf names
    :param newick: leave empty, this variable is used in recursion.
    :returns: tree in Newick format
    """
    
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parent_dist - node.dist, newick)
   
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parent_dist - node.dist, newick)
        else:
            newick = ");"
        newick = get_newick(node.get_left(), node.dist, leaf_names, newick=newick)
        newick = get_newick(node.get_right(), node.dist, leaf_names, newick=",%s" % (newick))
        newick = "(%s" % (newick)
        
        return newick
        
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
    
	# argument options
    
	parser = argparse.ArgumentParser(prog="partitioning_HC.py", formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent("""\
									###############################################################################             
									#                                                                             #
									#                             partitioning_HC.py                              #
									#                                                                             #
									###############################################################################  
									                            
									partitioning_HC.py obtains genetic clusters at any distance threshold(s) of an
									allele/SNP profile or pairwise distance matrix using hierarchical clustering 
									methods. If a profile is provided, pairwise hamming distances will be 
									calculated.
									
									Note: the profile matrix provided can only contain integers or IUPAC values.
									Alleles starting by "*" or "INF-" will be considered as a new. All the other 
									values will be considered as missing data.
									
									-----------------------------------------------------------------------------"""))
	
	group0 = parser.add_argument_group("Partitioning with Hierarchical Clustering", "Specifications to get partitions with HC methods")
	group0.add_argument("-d_mx", "--distance_matrix", dest="distance_matrix", required=False, type=str, help="[OPTIONAL] Input pairwise distance matrix")
	group0.add_argument("-a", "--allele-profile", dest="allele_profile", required=False, type=str, help="[OPTIONAL] Input allele profile matrix (can either be an allele matrix or a SNP matrix)")
	group0.add_argument("-l", "--loci", dest="loci", required=False, type=str, default = "none", help="[OPTIONAL] List of loci (e.g. cgMLST) that must be used for the clustering analysis. If \
					 	'--site-inclusion' argument > 0, this list of loci will be complemented with additional loci that fulfill the requirements specified in this argument.")
	group0.add_argument("-o", "--output", dest="out", required=True, type=str, help="[MANDATORY] Tag for output file name")
	group0.add_argument("--HC-threshold", dest="method_threshold", required=False, default="single", 
						help="[OPTIONAL] List of HC methods and thresholds to include in the analysis (comma-separated). To get clustering at all possible thresholds for a given method, just write \
						the method name (e.g. single). To get clustering at a specific threshold, indicate the threshold with a hyphen (e.g. single-10). To get clustering at a specific \
						range, indicate the range with a hyphen separating minimum and maximum (e.g. single-2-10). Note: Threshold values are inclusive, i.e. '--HC-threshold single-7' will consider \
						samples with <= 7 differences as belonging to the same cluster! Default: single (List of possible methods: single, complete, average, weighted, centroid, median, ward)")
	group0.add_argument("--pct-HC-threshold", dest="pct_HCmethod_threshold", required=False, default="none", help="[OPTIONAL] Similar to '--HC-threshold' but the partition threshold for cluster definition \
						is set as the proportion of differences to the final allelic schema size or number of informative positions, e.g. '--pct-HC-threshold single-0.005' corresponds to a \
						threshold of 5 allelic/SNP differences in a matrix with 1000 loci/sites under analysis. Ranges CANNOT be specified.")
	group0.add_argument("--site-inclusion", dest="samples_called", required=False, default = 0.0, help="[OPTIONAL: Useful to remove informative sites/loci with excess of missing data] Minimum \
						proportion of samples per site/loci without missing data (e.g. '--site-inclusion 1.0' will only keep loci/positions without missing data, i.e. a core alignment; \
						'--site-inclusion 0.0' will keep all loci/positions) NOTE: This argument works on profile/alignment loci/positions (i.e. columns)! [default: 0.0]")
	group0.add_argument("--loci-called", dest="loci_called", required=False, default = 0.0, help="[OPTIONAL] Minimum percentage of loci/positions called for allele/SNP matrices (e.g. \
						'--loci-called 0.95' will only keep in the profile matrix samples with > 95%% of alleles/positions, i.e. <= 5%% missing data). Applied after '--site-inclusion' argument! \
						[default: 0.0]")
	group0.add_argument("--missing-code", dest="missing_code", required=False, type=str, default = "0", help="[OPTIONAL] Code representing missing data. If different from '0' or 'N', please try \
		     			to avoid a IUPAC character (even in lower-case). [default: 0]")
	group0.add_argument("-m", "--metadata", dest="metadata", required=False, default="", type=str, help="[OPTIONAL] Metadata file in .tsv format to select the samples to use for clustering \
						according to the '--filter' argument")
	group0.add_argument("-f", "--filter", dest="filter_column", required=False, default="", help="[OPTIONAL] Filter for metadata columns to select the samples that must be used for HC \
						clustering. This must be specified within quotation marks in the following format 'column< >operation< >condition' (e.g. 'country == Portugal'). When \
						more than one condition is specified for a given column, they must be separated with commas (e.g 'country == Portugal,Spain,France'). When filters include more than one \
						column, they must be separated with semicolon (e.g. 'country == Portugal,Spain,France;date > 2018-01-01;date < 2022-01-01'). White spaces are important in this argument, \
						so, do not leave spaces before and after commas/semicolons.")
	group0.add_argument("-d", "--dist", dest="dist", required=False, default=1.0, type=float, help="Distance unit by which partition thresholds will be multiplied (example: if -d 10 and \
						--method-threshold single-2, the single linkage threshold will be set at 20).")
	
	
	args = parser.parse_args()


	# starting logs	----------

	log_name = args.out + ".log"
	log = open(log_name, "a+")
	
	print("\n-------------------- partitioning_HC.py --------------------\n")
	print("\n-------------------- partitioning_HC.py --------------------\n", file = log)
	print("version", version, "last updated on", last_updated, "\n")
	print("version", version, "last updated on", last_updated, "\n", file = log)
	print(" ".join(sys.argv), "\n")
	print(" ".join(sys.argv), "\n", file = log)
	
	missing_code = str(args.missing_code)
	if missing_code != "0":
		missing_need = True
	else:
		missing_need = False

	# processing allele profile ----------
	
	if args.allele_profile:
		print("Profile matrix provided... pairwise distance will be calculated!")
		print("Profile matrix provided... pairwise distance will be calculated!", file = log)
		
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
			missing_need = False
		
		# filtering allele matrix	----------
		
		if args.metadata != "" and args.filter_column != "":
			print("Filtering the allele matrix...")
			print("Filtering the allele matrix...", file = log)
			
			filters = args.filter_column
			mx = pandas.read_table(args.metadata, dtype = str)
			columns_names = [col.replace(" ", "_") for col in mx.columns]
			mx.replace(" ", "_", regex=True, inplace=True)
			mx.columns = columns_names
			sample_column = mx.columns[0]
			initial_samples = len(allele_mx[allele_mx.columns[0]].values.tolist())
			allele_mx = filter_mx(allele_mx, mx, filters, "allele", log)
			final_samples = len(allele_mx[allele_mx.columns[0]].values.tolist())
			print("\tFrom the " + str(initial_samples) + " samples, " + str(final_samples) + " were kept in the matrix...")
			print("\tFrom the " + str(initial_samples) + " samples, " + str(final_samples) + " were kept in the matrix...", file = log)
						
			if final_samples <= 1:
				print("\nCannot proceed because " + str(final_samples) + " samples were kept in the matrix!")
				print("\nCannot proceed because " + str(final_samples) + " samples were kept in the matrix!", file = log)
				sys.exit()
			allele_mx.to_csv(args.out + "_subset_matrix.tsv", sep = "\t", index = None)
	
	
		# cleaning allele matrix (columns)	----------
		pos_t0 = len(allele_mx.columns[1:])
		if args.loci != "none":
			loci2include = get_loci2use(args.loci,allele_mx)
			if float(args.samples_called) == 0.0:	
				print("Keeping the sites/loci present in the loci list.")
				print("Keeping the sites/loci present in the loci list.", file = log)			
				for col in allele_mx.columns[1:]:
					if col not in loci2include:
						allele_mx = allele_mx.drop(columns=col)
			else:
				print("Keeping the sites/loci present in the loci list and those with information in >= " + str(float(args.samples_called) * 100) + "%% of the samples...")
				print("Keeping the sites/loci present in the loci list and those with information in >= " + str(float(args.samples_called) * 100) + "%% of the samples...", file = log)
				for col in allele_mx.columns[1:]:
					if col not in loci2include:
						values = allele_mx[col].values.tolist()
						if (len(values)-values.count("0"))/len(values) < float(args.samples_called):
							allele_mx = allele_mx.drop(columns=col)
			pos_t1 = len(allele_mx.columns[1:])
			print("\tFrom the " + str(pos_t0) + " loci/positions, " + str(pos_t1) + " were kept in the matrix.")
			print("\tFrom the " + str(pos_t0) + " loci/positions, " + str(pos_t1) + " were kept in the matrix.", file = log)
		else:
			pos_t1 = len(allele_mx.columns[1:])
			if float(args.samples_called) > 0.0:
				print("Keeping the sites/loci with information in >= " + str(float(args.samples_called) * 100) + "%% of the samples...")
				print("Keeping the sites/loci with information in >= " + str(float(args.samples_called) * 100) + "%% of the samples...", file = log)
				for col in allele_mx.columns[1:]:
					values = allele_mx[col].values.tolist()
					if (len(values)-values.count("0"))/len(values) < float(args.samples_called):
						allele_mx = allele_mx.drop(columns=col)
				pos_t1 = len(allele_mx.columns[1:])
				print("\tFrom the " + str(pos_t0) + " loci/positions, " + str(pos_t1) + " were kept in the matrix.")
				print("\tFrom the " + str(pos_t0) + " loci/positions, " + str(pos_t1) + " were kept in the matrix.", file = log)
			
		if pos_t1 <= 1:
			print("\nCannot proceed because " + str(pos_t1) + " sites/loci were kept in the matrix!")
			print("\nCannot proceed because " + str(pos_t1) + " sites/loci were kept in the matrix!", file = log)
			sys.exit()
		with open(args.out + "_loci_used.txt", "w+") as loci_out:
			for locus in allele_mx.columns[1:]:
				print(locus, file = loci_out)
		allele_mx.to_csv(args.out + "_clean_missing_matrix.tsv", index = False, header=True, sep ="\t")
		
		# cleaning allele matrix (rows)	----------

		if float(args.loci_called) > 0.0 and float(args.samples_called) < 1.0:
			print("Cleaning the profile matrix using a threshold of >" + str(args.loci_called) + " alleles/positions called per sample...")
			print("Cleaning the profile matrix using a threshold of >" + str(args.loci_called) + " alleles/positions called per sample...", file = log)
			
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
			report_allele_df.to_csv(args.out + "_loci_report.tsv", index = False, header=True, sep ="\t")

			print("\tFrom the " + str(len(allele_mx[allele_mx.columns[0]].values.tolist())) + " samples, " + str(len(pass_samples)) + " were kept in the profile matrix.")
			print("\tFrom the " + str(len(allele_mx[allele_mx.columns[0]].values.tolist())) + " samples, " + str(len(pass_samples)) + " were kept in the profile matrix.", file = log)
			
			if len(pass_samples) <= 1:
				print("\nCannot proceed because " + str(len(pass_samples)) + " samples were kept in the matrix!")
				print("\nCannot proceed because " + str(len(pass_samples)) + " samples were kept in the matrix!", file = log)
				sys.exit()
			allele_mx = allele_mx[allele_mx[allele_mx.columns[0]].isin(pass_samples)]
			allele_mx.to_csv(args.out + "_flt_samples_matrix.tsv", index = False, header=True, sep ="\t")
			
			
		# getting distance matrix	----------
		
		print("Getting the pairwise distance matrix with cgmlst-dists (if your profile matrix is too big, this will be done in chunks of 2000 alleles/positions)...")
		print("Getting the pairwise distance matrix with cgmlst-dists (if your profile matrix is too big, this will be done in chunks of 2000 alleles/positions)...", file = log)
		
		# convert ATCG to integers
		alleles = conv_nucl(allele_mx, missing_code, missing_need)
		total_size = len(allele_mx.columns) - 1
		
		# divide a big dataframe into chunks
		start_chunk = 0
		chunk_size = 2000
		end_chunk = start_chunk + chunk_size
		df_counter = 0

		while end_chunk < total_size + chunk_size:
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
		dist = pandas.read_table(args.out + "_dist_hamming.tsv")
		
	elif args.distance_matrix:
		print("Distance matrix provided... pairwise distance will not be calculated!")
		print("Distance matrix provided... pairwise distance will not be calculated!", file = log)		
		dist = pandas.read_table(args.distance_matrix)
	
		# filtering the pairwise distance matrix	----------
		
		if args.metadata != "" and args.filter_column != "":
			print("Filtering the distance matrix...")
			print("Filtering the distance matrix...", file = log)
			
			filters = args.filter_column
			mx = pandas.read_table(args.metadata, dtype = str)
			columns_names = [col.replace(" ", "_") for col in mx.columns]
			mx.columns = columns_names
			initial_samples = len(dist[dist.columns[0]].values.tolist())
			sample_column = mx.columns[0]
			dist = filter_mx(dist, mx, filters, "dist", log)
			final_samples = len(dist[dist.columns[0]].values.tolist())
			print("\tFrom the " + str(initial_samples) + " samples, " + str(final_samples) + " were kept in the matrix...")
			print("\tFrom the " + str(initial_samples) + " samples, " + str(final_samples) + " were kept in the matrix...", file = log)
						
			if final_samples <= 1:
				print("\nCannot proceed because " + str(final_samples) + " samples were kept in the matrix!")
				print("\nCannot proceed because " + str(final_samples) + " samples were kept in the matrix!", file = log)
				sys.exit()
			dist.to_csv(args.out + "_flt_dist.tsv", sep = "\t", index = None)

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
	
	else:
		print("Could not find a profile or a distance matrix... One of them needs to be specified!!")
		print("Could not find a profile or a distance matrix... One of them needs to be specified!!", file = log)	
		sys.exit()
		
			
	# hierarchical clustering 	----------
	
	distance_mx, condensed_dist_mx, samples = dist_mx(dist, log)
	
	clustering = {"sequence": distance_mx.columns.tolist()}
	
	combinations2run = {}
	pct_correspondence = {}
	
	if args.pct_HCmethod_threshold != "none":
		print("\n\tCorrespondence between percentage and number of differences:")
		print("\n\tCorrespondence between percentage and number of differences:", file = log)
		print("\n\t#METHOD\tPERCENTAGE\tDIFFERENCES")
		print("\n\t#METHOD\tPERCENTAGE\tDIFFERENCES", file = log)
		for combination_pct in args.pct_HCmethod_threshold.split(","):
			method = combination_pct.split("-")[0]
			threshold_pct = combination_pct.split("-",1)[1]
			threshold = str(int(int(total_size) * float(threshold_pct)))
			info_run = threshold,"pct"
			
			if method not in combinations2run.keys():
				combinations2run[method] = []
			if info_run not in combinations2run[method]:
				combinations2run[method].append(info_run)

			if threshold not in pct_correspondence.keys():
				pct_correspondence[threshold] = []
			if str(threshold_pct) not in pct_correspondence[threshold]:
				pct_correspondence[threshold].append(str(threshold_pct))
				print("\t" + str(float(threshold_pct)*100) + "\t" + str(threshold))
				print("\t" + str(float(threshold_pct)*100) + "\t" + str(threshold), file = log)
		
	for combination in args.method_threshold.split(","):
		if "-" not in combination:
			threshold = "all"
			method = combination
			if method not in combinations2run.keys():
				combinations2run[method] = []
			info_run = threshold,"general"
			combinations2run[method].append(info_run)
		else:
			method = combination.split("-")[0]
			threshold = str(combination.split("-",1)[1])
			if method not in combinations2run.keys():
				combinations2run[method] = []
			info_run = threshold,"general"
			combinations2run[method].append(info_run)
	
	cluster_details = {}
	
	for method in combinations2run.keys():
		print("Hierarchical clustering with method: " + method + "...")
		print("Hierarchical clustering with method: " + method + "...", file = log)
		hc_matrix, max_dist = hcluster(condensed_dist_mx, method, log)
		
		# get newick
		
		print("\tGenerating newick file...")
		print("\tGenerating newick file...", file = log)
	
		tree = to_tree(hc_matrix, False)
		nw = get_newick(tree, tree.dist, samples)
		
		with open(args.out + "_" + method + "_HC.nwk", "w+") as newick_out:
			print(nw, file = newick_out)
		
		# partitioning
		
		print("\tDefining clusters...")
		print("\tDefining clusters...", file = log)
		
		
		for threshold,request in combinations2run[method]:
			if threshold == "all":
				print("\tCalculating clustering in range",str(0),str(max_dist),"with a distance of",str(args.dist))
				print("\tCalculating clustering in range",str(0),str(max_dist),"with a distance of",str(args.dist), file = log)
				for thr in range(0,int(max_dist) + 1):
					partition = method + "-" + str(thr) + "x" + str(args.dist)
					if partition not in cluster_details.keys():
						cluster_details[partition] = {}
					info_clusters = list(fcluster(hc_matrix, t = int(thr) * args.dist, criterion = "distance"))
					# change cluster name according to cluster size
					counter = {}
					singleton_counter = 0
					for cluster in set(info_clusters):
						counter[cluster] = info_clusters.count(cluster)
					for i in range(len(info_clusters)):
						if counter[info_clusters[i]] == 1:
							singleton_counter += 1
							info_clusters[i] = "singleton_" + str(singleton_counter)
							if info_clusters[i] not in cluster_details[partition].keys():
								cluster_details[partition][info_clusters[i]] = {}
								cluster_details[partition][info_clusters[i]][1] = []
							cluster_details[partition][info_clusters[i]][1].append(distance_mx.columns[i])
						else:
							cluster_size = counter[info_clusters[i]]
							info_clusters[i] = "cluster_" + str(info_clusters[i])
							if info_clusters[i] not in cluster_details[partition].keys():
								cluster_details[partition][info_clusters[i]] = {}
								cluster_details[partition][info_clusters[i]][cluster_size] = []
							cluster_details[partition][info_clusters[i]][cluster_size].append(distance_mx.columns[i])
					clustering[partition] = info_clusters
			else:
				if "-" in threshold:
					min_thr = int(threshold.split("-")[0])
					max_thr = int(threshold.split("-")[1]) + 1
					
					if min_thr < max_dist and max_thr > max_dist:
						max_thr = int(max_dist) + 1
					elif min_thr > max_dist:
						print("\tThe requested partition, " + str(min_thr) +  ", is higher than the maximum distance. Clustering will not be computed.")
						print("\tThe requested partition, " + str(min_thr) +  ", is higher than the maximum distance. Clustering will not be computed.", file = log)
						continue
					
					print("\tCalculating clustering in range",str(min_thr),str(max_thr),"with a distance of",str(args.dist))
					print("\tCalculating clustering in range",str(min_thr),str(max_thr),"with a distance of",str(args.dist), file = log)
					for thr in range(min_thr,max_thr):
						partition = method + "-" + str(thr) + "x" + str(args.dist)
						if partition not in cluster_details.keys():
							cluster_details[partition] = {}
						info_clusters = list(fcluster(hc_matrix, t = thr * args.dist, criterion = "distance"))
						# change cluster name according to cluster size
						counter = {}
						singleton_counter = 0
						for cluster in set(info_clusters):
							counter[cluster] = info_clusters.count(cluster)
						for i in range(len(info_clusters)):
							if counter[info_clusters[i]] == 1:
								singleton_counter += 1
								info_clusters[i] = "singleton_" + str(singleton_counter)
								if info_clusters[i] not in cluster_details[partition].keys():
									cluster_details[partition][info_clusters[i]] = {}
									cluster_details[partition][info_clusters[i]][1] = []
								cluster_details[partition][info_clusters[i]][1].append(distance_mx.columns[i])
							else:
								cluster_size = counter[info_clusters[i]]
								info_clusters[i] = "cluster_" + str(info_clusters[i])
								if info_clusters[i] not in cluster_details[partition].keys():
									cluster_details[partition][info_clusters[i]] = {}
									cluster_details[partition][info_clusters[i]][cluster_size] = []
								cluster_details[partition][info_clusters[i]][cluster_size].append(distance_mx.columns[i])
						clustering[partition] = info_clusters
				else:
					if request == "general":
						if int(threshold) > max_dist:
							print("\tThe requested partition, " + str(threshold) +  ", is higher than the maximum distance. Clustering will not be computed.")
							print("\tThe requested partition, " + str(threshold) +  ", is higher than the maximum distance. Clustering will not be computed.", file = log)
							continue
						partition = method + "-" + str(threshold) + "x" + str(args.dist)
						if partition not in cluster_details.keys():
							cluster_details[partition] = {}
						print("\tCalculating clustering for threshold",str(threshold),"with a distance of",str(args.dist))
						print("\tCalculating clustering for threshold",str(threshold),"with a distance of",str(args.dist), file = log)
						info_clusters = list(fcluster(hc_matrix, t = int(threshold) * args.dist, criterion = "distance"))
						# change cluster name according to cluster size
						counter = {}
						singleton_counter = 0
						for cluster in set(info_clusters):
							counter[cluster] = info_clusters.count(cluster)
						for i in range(len(info_clusters)):
							if counter[info_clusters[i]] == 1:
								singleton_counter += 1
								info_clusters[i] = "singleton_" + str(singleton_counter)
								if info_clusters[i] not in cluster_details[partition].keys():
									cluster_details[partition][info_clusters[i]] = {}
									cluster_details[partition][info_clusters[i]][1] = []
								cluster_details[partition][info_clusters[i]][1].append(distance_mx.columns[i])
							else:
								cluster_size = counter[info_clusters[i]]
								info_clusters[i] = "cluster_" + str(info_clusters[i])
								if info_clusters[i] not in cluster_details[partition].keys():
									cluster_details[partition][info_clusters[i]] = {}
									cluster_details[partition][info_clusters[i]][cluster_size] = []
								cluster_details[partition][info_clusters[i]][cluster_size].append(distance_mx.columns[i])
						clustering[partition] = info_clusters
					elif request == "pct":
						for percentage in pct_correspondence[threshold]:
							partition = method + "-" + str(percentage)
							if partition not in cluster_details.keys():
								cluster_details[partition] = {}
							print("\tCalculating clustering for threshold " + method + "-" + str(threshold) + ", which corresponds to the pct threshold of: " + str(percentage))
							print("\tCalculating clustering for threshold " + method + "-" + str(threshold) + ", which corresponds to the pct threshold of: " + str(percentage), file = log)
							info_clusters = list(fcluster(hc_matrix, t = int(threshold), criterion = "distance"))
							# change cluster name according to cluster size
							counter = {}
							singleton_counter = 0
							for cluster in set(info_clusters):
								counter[cluster] = info_clusters.count(cluster)
							for i in range(len(info_clusters)):
								if counter[info_clusters[i]] == 1:
									singleton_counter += 1
									info_clusters[i] = "singleton_" + str(singleton_counter)
									if info_clusters[i] not in cluster_details[partition].keys():
										cluster_details[partition][info_clusters[i]] = {}
										cluster_details[partition][info_clusters[i]][1] = []
									cluster_details[partition][info_clusters[i]][1].append(distance_mx.columns[i])
								else:
									cluster_size = counter[info_clusters[i]]
									info_clusters[i] = "cluster_" + str(info_clusters[i])
									if info_clusters[i] not in cluster_details[partition].keys():
										cluster_details[partition][info_clusters[i]] = {}
										cluster_details[partition][info_clusters[i]][cluster_size] = []
									cluster_details[partition][info_clusters[i]][cluster_size].append(distance_mx.columns[i])
							clustering[partition] = info_clusters
				

	# output partitions
	
	print("Creating sample partitions file...")
	print("Creating sample partitions file...", file = log)
	df_clustering = pandas.DataFrame(data = clustering)
	df_clustering.to_csv(args.out + "_partitions.tsv", sep = "\t", index = None)
    
    
    # output cluster composition
	
	print("Creating cluster composition file...")
	print("Creating cluster composition file...", file = log)
	cluster_composition = get_cluster_composition(args.out + "_clusterComposition.tsv", cluster_details)


	#output percentage correspondence

	if len(pct_correspondence.keys()) > 0:
		with open(args.out + "_parcentage2threshold.tsv", "w+") as out_pct:
			print("#percentage\tthreshold", file = out_pct)
			for threshold in pct_correspondence.keys():
				for percentage in pct_correspondence[threshold]:
					print(str(percentage) + "\t" + str(threshold), file = out_pct)

	print("\npartitioning_HC.py is done!")
	print("\npartitioning_HC.py is done!", file = log)

	log.close()

if __name__ == "__main__":
    main()
