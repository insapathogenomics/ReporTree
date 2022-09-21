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
import pandas
from datetime import date

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

partitioning_grapetree_script = os.path.realpath(__file__)
grapetree = partitioning_grapetree_script.rsplit("/", 1)[0] + "/GrapeTree/grapetree.py"

version = "1.0.0"
last_updated = "2022-08-26"

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
									
									
									NOTE: Do not forget to cite GrapeTree original authors.
									
									
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
group0.add_argument("-o", "--output", dest="out", required=True, type=str, help="[MANDATORY] Tag for output file name")
group0.add_argument("--site-inclusion", dest="samples_called", required=False, default = 0.0, help="[OPTIONAL: Useful to remove informative sites/loci with excess of missing data] Minimum \
					proportion of samples per site/loci without missing data (e.g. '--site-inclusion 1.0' will only keep loci/positions without missing data, i.e. a core alignment; \
					'--site-inclusion 0.0' will keep all loci/positions) NOTE: This argument works on profile/alignment loci/positions (i.e. columns)! [default: 1.0]. Code for missing data: 0.")
group0.add_argument("--loci-called", dest="loci_called", required=False, default = "", help="[OPTIONAL] Minimum percentage of loci/positions called for allele/SNP matrices (e.g. \
					'--loci-called 0.95' will only keep in the profile matrix samples with > 95%% of alleles/positions, i.e. <= 5%% missing data). Applied after '--site-inclusion' argument! \
					Code for missing data: 0.")
group0.add_argument("--method", dest="grapetree_method", default="MSTreeV2", help="\"MSTreeV2\" [DEFAULT]\n Alternative:\"MSTree (goeBURST)\"\n")
group0.add_argument("--missing", dest="handler", default=0, type=int, help="ONLY FOR MSTree. \n0: [DEFAULT] ignore missing data in pairwise comparison. \n1: remove column \
					with missing data. \n2: treat missing data as an allele. \n3: use absolute number of allelic differences.")
group0.add_argument("--wgMLST", dest="wgmlst", default=False, action="store_true", help="[EXPERIMENTAL] a better support of wgMLST schemes (check GrapeTree github for details).")
group0.add_argument("--n_proc",  dest="number_of_processes", type=int, default=5, help="Number of CPU processes in parallel use. [5]")
group0.add_argument("-thr", "--threshold", dest="threshold", default = "max", help="[OPTIONAL] Partition thresholds for clustering definition. Different thresholds can be comma-separated \
					(e.g. 5,8,16). Ranges can be specified with an hyphen (e.g. 5,8,10-20). If this option is not set, the script will perform clustering for all the values in the range 0 \
					to max. Note: Threshold values are inclusive, i.e. '-thr 7' will consider samples with <= 7 differences as belonging to the same cluster!")
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


# filtering allele matrix	----------

if args.metadata != "" and args.filter_column != "" and ".fasta" not in args.allele_profile and ".fa" not in args.allele_profile and ".fas" not in args.allele_profile:
	print("Filtering the allele matrix...")
	print("Filtering the allele matrix...", file = log)
	
	filters = args.filter_column
	mx = pandas.read_table(args.metadata, dtype = str)
	allele_mx = pandas.read_table(args.allele_profile, dtype = str)
	sample_column = mx.columns[0]
	
	if "date" in mx.columns and "iso_week" not in mx.columns:
		index_no = mx.columns.get_loc("date")
		mx["date"] = pandas.to_datetime(mx["date"], errors = "coerce")
		year = mx["date"].dt.isocalendar().year
		week = mx["date"].dt.isocalendar().week
		mx["iso_year"] = year.astype(str)
		mx["iso_week_nr"] = week.astype(str)
		mx["iso_week"] = year.astype(str).replace("<NA>", "-") + "W" + week.astype(str).replace("<NA>", "--")
		isoyear = mx.pop("iso_year")
		isoweek = mx.pop("iso_week_nr")
		isodate = mx.pop("iso_week")
		mx.insert(index_no + 1, "iso_year", isoyear)
		mx.insert(index_no + 2, "iso_week_nr", isoweek)
		mx.insert(index_no + 3, "iso_week", isodate)
	
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
	allele_mx_filtered.to_csv(args.out + "_subset_matrix.tsv", index = False, header=True, sep ="\t")
	allele_filename = args.out + "_subset_matrix.tsv"
	
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


# cleaning allele matrix (columns)

if args.samples_called != 0.0 and not os.path.exists(args.out + "_align_profile.tsv"):
	print("Keeping only sites/loci with information in >= " + str(float(args.samples_called) * 100) + "% of the samples...")
	print("Keeping only sites/loci with information in >= " + str(float(args.samples_called) * 100) + "% of the samples...", file = log)
	
	mx_allele = pandas.read_table(allele_filename, dtype = str)

	pos_t0 = len(mx_allele.columns[1:])
	for col in mx_allele.columns[1:]:
		values = mx_allele[col].values.tolist()
		if (len(values)-values.count("0"))/len(values) < float(args.samples_called):
			mx_allele = mx_allele.drop(columns=col)
		elif (len(values)-values.count(0))/len(values) < float(args.samples_called):
			mx_allele = mx_allele.drop(columns=col)
	mx_allele.to_csv(args.out + "_clean_missing_matrix.tsv", index = False, header=True, sep ="\t")
	allele_filename = args.out + "_clean_missing_matrix.tsv"
	pos_t1 = len(mx_allele.columns[1:])
	print("\tFrom the " + str(pos_t0) + " loci/positions, " + str(pos_t1) + " were kept in the matrix.")
	print("\tFrom the " + str(pos_t0) + " loci/positions, " + str(pos_t1) + " were kept in the matrix.", file = log)
		

# cleaning allele matrix (rows)	----------

if args.loci_called != "":
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
	flt_report = report_allele_df[report_allele_df["pct_called"] > float(args.loci_called)]
	pass_samples = flt_report["samples"].values.tolist()
	
	print("\tFrom the " + str(len(allele_mx[allele_mx.columns[0]].values.tolist())) + " samples, " + str(len(pass_samples)) + " were kept in the profile matrix.")
	print("\tFrom the " + str(len(allele_mx[allele_mx.columns[0]].values.tolist())) + " samples, " + str(len(pass_samples)) + " were kept in the profile matrix.", file = log)
	
	allele_mx = allele_mx[allele_mx[allele_mx.columns[0]].isin(pass_samples)]
	allele_mx.to_csv(args.out + "_flt_samples_matrix.tsv", index = False, header=True, sep ="\t")
	allele_filename = args.out + "_flt_samples_matrix.tsv"
	report_allele_df.to_csv(args.out + "_loci_report.tsv", index = False, header=True, sep ="\t")
	
	
# preparing allele matrix for grapetree	----------

if args.matrix4grapetree:
	mx_allele = pandas.read_table(allele_filename, dtype = str)
	first_col = str(mx_allele.columns[0])
	if first_col[0] != "#":
		mx_allele = mx_allele.rename(columns={first_col: "#" + first_col})
		mx_allele.to_csv(args.out + "_alleles_4_grapetree.tsv", index = False, header=True, sep ="\t")
		
		
# running grapetree	----------

print("Running GrapeTree...")
print("Running GrapeTree...", file = log)

if args.wgmlst == True:
	print("python " + grapetree + " -p " + allele_filename + " -m " + args.grapetree_method + " -o " + args.out + " --missing " + str(args.handler) + " --n_proc " + str(args.number_of_processes) + " --wgMLST")
	print("python " + grapetree + " -p " + allele_filename + " -m " + args.grapetree_method + " -o " + args.out + " --missing " + str(args.handler) + " --n_proc " + str(args.number_of_processes) + " --wgMLST", file = log)
	os.system("python " + grapetree + " -p " + allele_filename + " -m " + args.grapetree_method + " -o " + args.out + " --missing " + str(args.handler) + " --n_proc " + str(args.number_of_processes) + " --wgMLST")
else:
	print("python " + grapetree + " -p " + allele_filename + " -m " + args.grapetree_method + " -o " + args.out + " --missing " + str(args.handler) + " --n_proc " + str(args.number_of_processes))
	print("python " + grapetree + " -p " + allele_filename + " -m " + args.grapetree_method + " -o " + args.out + " --missing " + str(args.handler) + " --n_proc " + str(args.number_of_processes), file = log)
	os.system("python " + grapetree + " -p " + allele_filename + " -m " + args.grapetree_method + " -o " + args.out + " --missing " + str(args.handler) + " --n_proc " + str(args.number_of_processes))


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

for min_r,max_r in range_def:
	print("\tCalculating clustering in range",min_r,max_r,"with a distance of",distance_magnitude)
	print("\tCalculating clustering in range",min_r,max_r,"with a distance of",distance_magnitude, file = log)
	for partition_value in range(min_r,max_r):
		partition = partition_value * distance_magnitude
		if partition <= max_distance:
			order_partitions.append("MST-" + str(partition_value) + "x" + str(distance_magnitude))
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

		info["MST-" + str(partition_value) + "x" + str(distance_magnitude)] = clusters


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
	for min_r,max_r in range_def:
		for partition_value in range(min_r,max_r):
			partition = partition_value * distance_magnitude
			if "MST-" + str(partition_value) + "x" + str(distance_magnitude) not in typing.keys():
				typing["MST-" + str(partition_value) + "x" + str(distance_magnitude)] = []
			typing["MST-" + str(partition_value) + "x" + str(distance_magnitude)].append(info_sample[sample]["MST-" + str(partition_value) + "x" + str(distance_magnitude)])

matrix = pandas.DataFrame(data = typing, columns = order_partitions)
matrix.to_csv(args.out + "_partitions.tsv", index = False, header=True, sep ="\t")


print("\npartitioning_grapetree.py is done!")
print("\npartitioning_grapetree.py is done!", file = log)

log.close()
