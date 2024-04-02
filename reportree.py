#!/usr/bin/env	python3

"""
Obtain clustering information with ReporTree

By Veronica Mixao
@INSA
"""

import os
import sys
import argparse
import textwrap
import datetime as datetime
from datetime import date
import pandas

version = "2.4.1"
last_updated = "2024-04-02"

reportree_script = os.path.realpath(__file__)
reportree_path = reportree_script.rsplit("/", 1)[0]
python = sys.executable


# command line functions	----------

def run_stability(args, log_name, partitions_status):
	""" run comparing_partitions_v2.py """
	
	extra_stability = ""
	if args.keep_redundants:
		extra_stability += " --keep-redundants "
	if partitions_status == "new":
		partitions = args.output + "_partitions.tsv"
	elif partitions_status == "old":
		partitions = args.partitions
	cmd = python + " " + reportree_path + "/scripts/ComparingPartitions/comparing_partitions_v2.py -i1 " + partitions + " -t " + args.output + " -o1 " +  str(args.order) + " -a stability -n \
	" + str(args.n_obs) + " -thr " + str(args.AdjustedWallace) + " -log " + log_name + " " + extra_stability
	returned_value = os.system(cmd)

	return returned_value

def run_metadata_report(args, partitions2report_final, partitions_status):
	""" run metadata_report.py """

	extra_metadata = ""
	if args.mx_transpose:
		extra_metadata += " --mx-transpose "
	if args.pivot:
		extra_metadata += " --pivot "
	if partitions_status == "new":
		partitions = args.output + "_partitions.tsv"
		no_partitions = False
	elif partitions_status == "old":
		partitions = args.partitions
		no_partitions = False
	else:
		no_partitions = True
	if no_partitions:
		cmd = python + " " + reportree_path + "/scripts/metadata_report.py -m " + args.metadata + " -o " + args.output + " --columns_summary_report \
		" + args.columns_summary_report + " --partitions2report " + partitions2report_final + " --metadata2report " + args.metadata2report + " -f \"" + args.filter_column + "\" \
		--frequency-matrix \'" + args.frequency_matrix + "\' --count-matrix \'" + args.count_matrix + "\' " + extra_metadata	
	else:
		cmd = python + " " + reportree_path + "/scripts/metadata_report.py -m " + args.metadata + " -p " + partitions + " -o " + args.output + " --columns_summary_report \
		" + args.columns_summary_report + " --partitions2report " + partitions2report_final + " --metadata2report " + args.metadata2report + " -f \"" + args.filter_column + "\" \
		--frequency-matrix \'" + args.frequency_matrix + "\' --count-matrix \'" + args.count_matrix + "\' " + extra_metadata
	returned_value = os.system(cmd)

	return returned_value

def run_treecluster(args):
	""" run partitioning_treecluster.py """
	
	extra_treecluster = ""
	if not args.root_dist_by_node:
		extra_treecluster += " --root-dist-by-node "
	if args.support != "-inf":
		extra = " --support " + str(args.support)
		extra_treecluster += extra
	if args.subset == True:
		extra = " -m " + args.metadata + " -f \"" + args.filter_column + "\""
		extra_treecluster += extra
	cmd = python + " " + reportree_path + "/scripts/partitioning_treecluster.py -t " + args.tree + " -o " + args.output + " -d " + str(args.dist) + " --method-threshold \
	" + args.method_threshold + " --root " + args.root + " " + extra_treecluster
	returned_value = os.system(cmd)
	
	return returned_value

def run_alignment_processing(args):
	""" run alignment_processing.py """
	
	extras_alignment = ""
	if args.subset:
		if "\"" in args.filter_column:
			args.filter_column = str(args.filter_column).split("\"")[1]
		extra = " -m " + args.metadata + " -f \"" + args.filter_column + "\""
		extras_alignment += extra
	if args.remove_ref:
		extras_alignment += " --remove-reference "
	if args.use_ref:
		extras_alignment += " --use-reference-coords "
	if args.remove_ref or args.use_ref:
		extras_alignment += " -r " + args.reference
	if args.missing_code == "0":
		missing_code = "N"
	else:
		missing_code = args.missing_code 
	cmd = python + " " + reportree_path + "/scripts/alignment_processing.py -align " + args.alignment + " -o " + args.output + " --sample-ATCG-content " + str(args.ATCG_content) + " \
	-r " + args.reference + " --site-ATCG-content " + str(args.N_content) + " --missing-code " + str(missing_code) + " --get-position-correspondence " + args.pos_corr + " --position-list \
	" + args.pos_list + " " + extras_alignment
	returned_value = os.system(cmd)

	return returned_value

def run_vcf2mst(args, intype):
	""" run vcf2mst.pl """

	if intype == "vcf":
		cmd = "perl " + reportree_path + "/scripts/vcf2mst/vcf2mst.pl " + args.vcf + " " + args.output + "_profile.tsv vcf -out profile"
	elif intype == "var":
		cmd = "perl " + reportree_path + "/scripts/vcf2mst/vcf2mst.pl " + args.variants + " " + args.output + "_profile.tsv tsv -out profile -tsv-sample-pos 0 \
		-tsv-mutationslist-find pos -tsv-mutationslist-pos 1 -tsv-mutation-pos-regexp '\w(\d+)\w'"
	returned_value = os.system(cmd)

	return returned_value

def run_partitioning_grapetree(args, input_align, profile):
	""" run partitioning_grapetree.py """

	extras_grapetree = ""
	if args.matrix4grapetree == True:
		extras_grapetree += " --matrix-4-grapetree "
	if args.wgmlst == True:
		extras_grapetree += " --wgMLST "
	if args.subset == True:
		if not input_align:
			extra = " -m " + args.metadata + " -f \"" + args.filter_column + "\""
			extras_grapetree += extra
	if input_align:
		site_inclusion = "0"
		loci_called = "0"
		missing_code = "0"
	else:
		site_inclusion = str(args.N_content)
		loci_called = str(args.loci_called)
		missing_code = str(args.missing_code)
	cmd = python + " " + reportree_path + "/scripts/partitioning_grapetree.py -a " + profile + " -o " + args.output + " --method " + args.grapetree_method + " --missing \
	" + str(args.handler) + " --n_proc " + str(args.number_of_processes) + " -thr " + str(args.threshold) + " -d " + str(args.dist) + " --site-inclusion " + site_inclusion + " \
	-pct_thr " + str(args.pct_threshold) + " --loci-called " + loci_called + " --missing-code " + missing_code + " " + extras_grapetree
	returned_value = os.system(cmd)

	return returned_value

def run_partitioning_HC(args, input_align, distance_matrix_input, profile):
	""" run partitioning_HC.py """
	
	extras_hc = ""
	if args.subset == True:
		if not input_align:
			extra = " -m " + args.metadata + " -f \"" + args.filter_column + "\""
			extras_hc += extra		
	if not distance_matrix_input:
		if input_align:
			site_inclusion = "0"
			loci_called = "0"
			missing_code = "0"
		else:
			site_inclusion = str(args.N_content)
			loci_called = str(args.loci_called)
			missing_code = str(args.missing_code)
		cmd = python + " " + reportree_path + "/scripts/partitioning_HC.py -a " + profile + " -o " + args.output + " --HC-threshold " + args.HCmethod_threshold + " \
		-d " + str(args.dist) + " --site-inclusion " + site_inclusion + " --loci-called " + loci_called + " --missing-code " + missing_code + " --pct-HC-threshold \
		" + str(args.pct_HCmethod_threshold) + " " + extras_hc
	else:
		cmd = python + " " + reportree_path + "/scripts/partitioning_HC.py -d_mx " + profile + " -o " + args.output + " --HC-threshold " + args.HCmethod_threshold + " \
		-d " + str(args.dist) + " --pct-HC-threshold " + str(args.pct_HCmethod_threshold) + " " + extras_hc
	returned_value = os.system(cmd)
	
	return returned_value


# basic functions	----------

def print_log(message, log):
	""" print messages in the terminal and in the log file """

	print(message)
	print(message, file = log)
	
def infer_partitions_from_list(analysis, partitions, metadata, thresholds, output, distance, log):
	""" infer which partition names should be provided for a given analysis from a list of input thresholds 
	input: list of thresholds provided in the command lines
	output: list of partition names """

	distance_code = "x" + str(distance)
	partitions_mx = pandas.read_table(partitions)

	if thresholds == "all":
		partitions2report = "all"
		partitions2report_lst = ["all"]
	else:
		partitions2report_lst = []
		for part in thresholds.split(","):
			if part in metadata:
				partitions2report_lst.append(part)
			elif part == "stability_regions": # get stability regions
				with open(output + "_stableRegions.tsv", "r") as stable_blocks_open:
					stable_blocks = stable_blocks_open.readlines()
					for line in stable_blocks:
						if "#" not in line:
							l = line.split("\t")
							first_partition = l[1].split("->")[1]
							if first_partition not in partitions2report_lst:
								partitions2report_lst.append(first_partition)
			elif "nomenclature_code" in part: # this is the nomenclature code
				partitions2report_lst.append(part)
			else:
				if analysis == "grapetree":
					method_code = "MST-"
					if "-" in part: # range specified
						min_threshold = part.split("-")[0]
						max_threshold = part.split("-")[1]
						if str(min_threshold).isnumeric() and str(max_threshold).isnumeric():
							for threshold in range(int(min_threshold),int(max_threshold) + 1):
								info = method_code + str(threshold) + distance_code
								partitions2report_lst.append(info)
						else:
							print_log("\tOnly integer values are allowed for range specification. We are sorry but your request " + str(part) + " will be ignored.", log)
					else: # partition corresponds to a single value
						if "." not in part:
							if str(part).isnumeric():
								info = method_code + str(part) + distance_code
								partitions2report_lst.append(info)
						elif "." in part: # it is a percentage
							for col in partitions_mx.columns:
								if method_code in col and str(part) in col:
									partitions2report_lst.append(col)
						else:
							print_log("\tWe are sorry but we cannot solve your request " + str(part) + ", so, it will be ignored.", log)
				elif analysis == "treecluster" or analysis == "HC":
					if "-" in part:
						method_code = part.split("-")[0] + "-"
						threshold_info = part.split("-", 1)[1]
						if "-" in threshold_info: # range specified
							min_threshold = threshold_info.split("-")[0]
							max_threshold = threshold_info.split("-")[1]
							if str(min_threshold).isnumeric() and str(max_threshold).isnumeric():
								for threshold in range(int(min_threshold),int(max_threshold) + 1):
									info = method_code + str(threshold) + distance_code
									partitions2report_lst.append(info)
							else:
								print_log("\tOnly integer values are allowed for range specification. We are sorry but your request " + str(part) + " will be ignored.", log)
						else: # partition corresponds to a single value
							if str(threshold_info).isnumeric(): # it is an integer
								info = method_code + str(threshold_info) + distance_code
								partitions2report_lst.append(info)
							elif "." in threshold_info: # it is a percentage
								for col in partitions_mx.columns:
									if str(part) == col:
										partitions2report_lst.append(col)
							else:
								print_log("\tWe are sorry but we cannot solve your request " + str(part) + ", so, it will be ignored.", log)
				elif analysis == "other":
					partitions2report_lst.append(part)
	
		partitions2report = ",".join(partitions2report_lst)
	
	return partitions2report, partitions2report_lst

def all_partitions_available(method_threshold, include_node):
	""" determines if all partitions are available for any method,
	and, if so, which methods 
	
	returns:
	list of methods with all partitions
	"""
	
	methods =[]
	 
	for mt in method_threshold.split(","):
		if "-" not in mt:
			methods.append(mt)
	
	if include_node == False:
		methods.append("node")
	
	return methods
			
def filter_partitions_table(method, table, output):
	""" only keep the partitions table of the method in analysis """
	
	mx = pandas.read_table(table)
	
	suitable_columns = [mx.columns[0]]
	for col in mx.columns.tolist():
		if col.split("-")[0] == method:
			suitable_columns.append(col)
	
	final_df = mx.filter(items=suitable_columns)
	
	with open(output + "_tmp.tsv", "w") as method_file:
		final_df.to_csv(output + "_tmp.tsv", index = False, header=True, sep ="\t")
		
def col_list(args):
	""" output the list of columns 
	input: argument list
	returns: list with the column names
	"""
	
	cols_output = []
	methods_output = []
	metadata_mx = pandas.read_table(args.metadata, dtype = str)
	columns_names = [col.replace(" ", "_") for col in metadata_mx.columns]
	metadata_mx.columns = columns_names
	
	for col in metadata_mx.columns:
		if col != "date":
			n_col = "n_" + col
			cols_output.append(col)
			cols_output.append(n_col)
		else:
			cols_output.append("first_seq_date")
			cols_output.append("last_seq_date")
			cols_output.append("timespan_days")
	
	if args.partitions != "":
		partitions_mx = pandas.read_table(args.partitions)
		columns_names = [col.replace(" ", "_") for col in partitions_mx.columns]
		partitions_mx.columns = columns_names
		for col in partitions_mx.columns:
			if col not in cols_output:
				n_col = "n_" + col
				cols_output.append(col)
				cols_output.append(n_col)
	
	elif args.tree != "":
		method_threshold = args.method_threshold
		if "," in method_threshold:
			for mt in method_threshold.split(","):
				method = mt.split("-")[0]
				if method not in methods_output:
					methods_output.append(method)
		else:
			method = method_threshold.split("-")[0]
			methods_output.append(method)
			
	elif args.distance_matrix != "" or args.analysis == "HC":
		HC_method_threshold = args.HCmethod_threshold
		if "," in HC_method_threshold:
			for mt in HC_method_threshold.split(","):
				method = mt.split("-")[0]
				if method not in methods_output:
					methods_output.append(method)
		else:
			method = HC_method_threshold.split("-")[0]
			methods_output.append(method)	
	
	elif args.analysis == "grapetree":
		methods_output.append("MST")
	
	return cols_output,methods_output

def info_samples_interest(samples, matrix, partitions, out):
	""" filter partitions summary """
	
	if "," in samples:
		samples_of_interest = samples.split(",")
	elif ".tsv" in samples or ".txt" in samples or ".csv" in samples:
		samples_of_interest = []
		with open(samples, "r") as input_file:
			infile = input_file.readlines()
			for line in infile:
				l = line.split("\n")[0]
				if "," in l:
					info = l.split(",")
					for s in info:
						samples_of_interest.append(s)
				elif "\t" in l:
					info = l.split("\t")[0]
					samples_of_interest.append(info)
				else:
					samples_of_interest.append(l)
	else:
		samples_of_interest = []
		samples_of_interest.append(samples)
				
	found = set()
	partitions_mx = pandas.read_table(partitions)
	singletons = []
	do_not_exist = []

	with open(matrix, "r") as mx:
		with open(out + "_SAMPLES_OF_INTEREST_partitions_summary.tsv", "w+") as outfile:
			i = 0
			for line in mx.readlines():
				lin = line.split("\n")
				l = lin[0].split("\t")
				
				if i == 0:
					print("SAMPLE_OF_INTEREST\t" + lin[0], file = outfile)
					if "nomenclature_change" in line:
						nomenclature = True
					else:
						nomenclature = False
				else:
					if nomenclature:
						samples_cluster = l[5].split(",")
					else:
						samples_cluster = l[3].split(",")
					samples_observed = []
					for sample in samples_of_interest:
						if sample in samples_cluster:
							samples_observed.append(sample)
							found.add(sample)
					if len(samples_observed) >= 1:
						print(",".join(samples_observed) + "\t" + lin[0], file = outfile)
				i += 1
			
			if len(found) != len(samples_of_interest):
				not_found = set(samples_of_interest) - found
				samples_analyzed = list(partitions_mx[partitions_mx.columns[0]])
				
				for sample in not_found:
					if sample in samples_analyzed:
						singletons.append(sample)
					else:
						do_not_exist.append(sample)
				
				if len(singletons) > 0:
					print("*Sample(s) " + ",".join(singletons) + " are singletons at all thresholds used!", file = outfile)
				if len(do_not_exist) > 0:
					print("**Sample(s) " + ",".join(do_not_exist) + " were not found in the partitions table! Please check if they were provided in the input file.", file = outfile)

	return samples_of_interest, singletons, do_not_exist

def interest2metadata(tag, samples_of_interest):
	""" update metadata_w_partitions with information about samples of interest from the current run 
	input: list of samples and metadata table
	output: metadata table """
	
	mx = pandas.read_table(tag + "_metadata_w_partitions.tsv")
	samples = mx[mx.columns[0]].values.tolist()
	category = []
	for s in samples:
		if s in samples_of_interest:
			category.append("sample of interest")
		else:
			category.append("")
	if "category" in mx.columns:
		index_no = mx.columns.get_loc("category")
		mx["category_original"] = mx["category"]
		category_original = mx.pop("category_original")
		mx.insert(index_no, "category_original", category_original)
		mx = mx.drop(["category"], axis=1)
	
	i = 0
	for col in mx.columns:
		if "single-" in col or "MST-" in col:
			break
		i += 1			
	mx.insert(i, "category", category, True)
	mx.to_csv(tag + "_metadata_w_partitions.tsv", index = False, header=True, sep = "\t")
	metadata = tag + "_metadata_w_partitions.tsv"
	
	return metadata
	
def get_clusters_interest(samples, matrix, metadata, day, partitions):
	""" get clusters of interest 
	input: list of samples and partitions table
	output: list of clusters of interest """

	partitions_mx = pandas.read_table(matrix)
	metadata_mx = pandas.read_table(metadata, dtype = str)
	nomenclature_code = "nomenclature_code_" + day
	clusters_of_interest = {}
	not_in_partitions = []

	for sample in samples:
		flt_mx = partitions_mx[partitions_mx[partitions_mx.columns[0]] == sample]
		for p in partitions:
			if p in partitions_mx.columns:
				if len(flt_mx[p].values.tolist()) > 0:
					cluster = flt_mx[p].values.tolist()[0]
					cluster_length = len(partitions_mx[partitions_mx[p] == cluster].values.tolist())
					if cluster_length > 1:
						if p not in clusters_of_interest.keys():
							clusters_of_interest[p] = []
						if cluster not in clusters_of_interest[p]:
							clusters_of_interest[p].append(cluster)
			elif p == nomenclature_code and p in metadata_mx.columns:
				flt_metadata = metadata_mx[metadata_mx[metadata_mx.columns[0]] == sample]
				if len(flt_metadata[p].values.tolist()) > 0:
					cluster = flt_metadata[p].values.tolist()[0]
					cluster_length = len(metadata_mx[metadata_mx[p] == cluster].values.tolist())
					if cluster_length > 1:
						if p not in clusters_of_interest.keys():
							clusters_of_interest[p] = []
						if cluster not in clusters_of_interest[p]:
							clusters_of_interest[p].append(cluster)
			else:
				not_in_partitions.append(p)			

	return clusters_of_interest, not_in_partitions


def loci_called2metadata(metadata_original, out, thr, analysis):
	""" adds column with percentage of called loci
	to metadata table """
	
	metadata = pandas.read_table(metadata_original, dtype = str)
	columns_names = [col.replace(" ", "_") for col in metadata.columns]
	metadata.columns = columns_names
	loci = pandas.read_table(out + "_" + analysis + "_report.tsv")
	a = metadata.set_index(metadata.columns[0], drop = True)
	b = loci.set_index(loci.columns[0], drop = True)
	
	b["QUAL_called"] = b["pct_called"]
	b["QUAL_called"].where(b["QUAL_called"] <= float(thr), "PASS", inplace = True)
	b["QUAL_called"].where(b["QUAL_called"].astype(str) == "PASS", "excluded", inplace = True)

	c = pandas.concat([a, b["pct_called"], b["QUAL_called"]], axis=1)
	c = c.reset_index(drop = False)
	c.rename(columns={c.columns[0]: metadata.columns[0]}, inplace=True)
	
	complete_metadata = pandas.DataFrame(data = c)
	complete_metadata.to_csv(metadata_original, index = False, header=True, sep = "\t")

def get_closest_samples(sample, n, hamming):
	""" get the list of n closest samples for subtree 
	input: sample, n and hamming matrix 
	output: list of samples """

	d_mx = pandas.read_table(hamming)
	if sample not in d_mx[d_mx.columns[0]].values.tolist():
		code = "not_in_hamming"
		closest_samples = ""
	else:
		d_mx.sort_values(by=[sample], inplace = True)
		distances = d_mx[sample].values.tolist()
		samples = d_mx[d_mx.columns[0]].values.tolist()
		if int(n) - 1 < len(distances):
			max_dist = distances[int(n)]
			closest_samples = []
			i = 0
			for d in distances:
				if d <= max_dist:
					closest_samples.append(samples[i])
					i += 1
			code = "run"
		else:
			closest_samples = samples
			code = "all_samples"

	return closest_samples, code

def metadata4subset(mx, partitions_mx, tag_subset, filter_subset):
	""" generate a new metadata table only with the subset of samples
	input: big metadata, tag_subset, filter_subset
	output: 'subset.tsv' """

	mx = pandas.read_table(mx)
	
	filter_subset = filter_subset.split(" ")
	col = filter_subset[0]
	val = filter_subset[2].split(",")
	if col not in mx.columns:
		partitions_mx = pandas.read_table(partitions_mx)
		partitions_mx.set_index(partitions_mx.columns[0], drop = False, inplace = True)
		mx.set_index(mx.columns[0], drop = False, inplace = True)
		if col in partitions_mx.columns:
			mx = pandas.concat([mx, partitions_mx[col]], axis=1)
	new_mx = mx.loc[mx[col].isin(val)]
	mx_name = tag_subset + "_metadata.tsv"
	new_mx.to_csv(mx_name, index = False, header=True, sep ="\t")

	return mx_name	 

def get_method_threshold(method_threshold):
	""" infer the method threshold for subset when running HC """
	
	final_methods = []
	method_threshold = method_threshold.split(",")
	for m_thr in method_threshold:
		if "-" not in m_thr:
			method = m_thr
		else:
			method = m_thr.split("-")[0]
		if method not in final_methods:
			final_methods.append(method)
	
	return ",".join(final_methods)


# functions cluster nomenclature	----------

def filter_db(mx_original, mx_new):
    """ Keep only the samples of the run 
    input: dataframe partitions, dataframe nomenclature
    output: filtered dataframe nomenclature """

    samples = mx_original[mx_original.columns[0]].values.tolist()
    mx_new = mx_new[mx_new[mx_new.columns[0]].isin(samples)]
    old_samples = mx_new[mx_new.columns[0]].values.tolist()

    return mx_new, old_samples

def get_new_db(mx_new):
    """ Generate a new dictionary with the new cluster names 
    input: pandas dataframe
    output: dictionary """
    
    db_samples = {}
    
    for classification in mx_new.columns[1:]:
        clusters = mx_new[classification].values.tolist()
        db_samples[classification] = {}
        for cluster in set(clusters):
            cluster_mx = mx_new.loc[mx_new[classification] == cluster]
            samples = ",".join(cluster_mx[cluster_mx.columns[0]].values.tolist())
            db_samples[classification][samples] = str(cluster)
    
    return db_samples

def get_last_counters(clusters):
    """ obtain the last number used for clusters and singletons by ReporTree """

    cluster_val = []
    singleton_val = []
    clusters = clusters.values.tolist()

    for observation in set(clusters):
        if "cluster_" in observation:
            info = observation.split("cluster_")[1]
            if "." in info:
                n = info.split(".")[0]
            else:
                n = info
            if n.isdigit():
                cluster_val.append(n)
        elif "singleton_" in observation:
            info = observation.split("singleton_")[1]
            n = info
            if n.isdigit():
                singleton_val.append(n)
    
    if len(cluster_val) >= 1:
        max_cluster = max(sorted(map(int,cluster_val)))
    else:
        max_cluster = 0
    if len(singleton_val) >= 1:
        max_singleton = max(sorted(map(int,singleton_val)))
    else:
        max_singleton = 0

    return int(max_cluster), int(max_singleton)

def determine_new_cluster_names(previous_samples, db_samples, old_partition, mx_new, partition, mx_original, max_cluster, max_singleton):
    """ Determine the cluster name conversion
    input: dictionary with wanted names and dataframe
    output: dictionary with conversion """
    
    conversion = {} # cluster_names_db[cluster partitions] = cluster nomenclature
    mx = mx_original.set_index(mx_original.columns[0], drop = False)
    clusters = list(mx[partition].unique())
    changes = {}
    split_clusters = {}

    for cluster in clusters:
        merged_clusters = {}
        cluster_mx = mx[mx[partition] == cluster]
        cluster_samples = cluster_mx[cluster_mx.columns[0]].values.tolist()
        new_samples = set(cluster_samples) - set(previous_samples)
        if len(new_samples) == len(cluster_samples): # all cluster samples are new
            if len(new_samples) == 1:
                max_singleton += 1
                new_cluster_name = "singleton_" + str(max_singleton)
            else:
                max_cluster += 1
                new_cluster_name = "cluster_" + str(max_cluster)
            conversion[cluster] = new_cluster_name
            change_info = "new","-","new",new_cluster_name,str(len(cluster_samples)),str(len(new_samples)),",".join(list(new_samples))
            if "new" not in changes.keys():
                changes["new"] = []
            changes["new"].append(change_info)
        elif len(new_samples) == 0: # cluster does not have new samples
            for sample_list in db_samples.keys():
                old_samples = sample_list.split(",")
                old_cluster = str(db_samples[sample_list])
                cluster_intersection = set(cluster_samples) & set(old_samples)
                if len(cluster_intersection) > 0: # they share something
                    exclusive_old = set(old_samples) - set(cluster_samples)
                    exclusive_new = set(cluster_samples) - set(old_samples)
                    if len(exclusive_old) == 0 and len(exclusive_new) == 0: # cluster remains exactly the same
                        conversion[cluster] = old_cluster
                        change_info = old_cluster,str(len(old_samples)),"kept (none)",old_cluster,str(len(cluster_samples)),str(len(new_samples)),",".join(list(new_samples))
                        if old_cluster not in changes.keys():
                            changes[old_cluster] = []
                        changes[old_cluster].append(change_info)
                    elif len(exclusive_old) > 0 and len(exclusive_new) == 0: # cluster was split and did not merge
                        if old_cluster not in split_clusters.keys():
                            split_clusters[old_cluster] = 0
                        if old_cluster not in changes.keys():
                            changes[old_cluster] = []
                        if len(cluster_samples) == 1: # it is a singleton now
                            max_singleton += 1
                            new_cluster_name = "singleton_" + str(max_singleton)
                            conversion[cluster] = new_cluster_name
                            change_info = old_cluster,str(len(old_samples)),"new (split)",new_cluster_name,str(len(cluster_samples)),str(len(new_samples)),",".join(list(new_samples))
                            changes[old_cluster].append(change_info)
                        else: # it is a cluster
                            split_clusters[old_cluster] += 1
                            new_cluster_name = old_cluster + "." + str(split_clusters[old_cluster])
                            conversion[cluster] = new_cluster_name
                            change_info = old_cluster,str(len(old_samples)),"new (split)",new_cluster_name,str(len(cluster_samples)),str(len(new_samples)),",".join(list(new_samples))
                            changes[old_cluster].append(change_info)
                    elif len(exclusive_old) > 0 and len(exclusive_new) > 0: # cluster was split and did merge
                        cluster_old_mx = mx_new[mx_new[mx_new.columns[0]].isin(list(cluster_samples))]
                        all_clusters = sorted(cluster_old_mx[old_partition].unique())
                        all_clusters = "_".join(all_clusters)
                        if all_clusters not in merged_clusters.keys():
                            max_cluster += 1
                            new_cluster_name = "cluster_" + str(max_cluster)
                            merged_clusters[all_clusters] = new_cluster_name
                        else:
                            new_cluster_name = merged_clusters[all_clusters]
                        conversion[cluster] = new_cluster_name
                        change_info = all_clusters,"multiple clusters","new (split_merge)",new_cluster_name,str(len(cluster_samples)),str(len(new_samples)),",".join(list(new_samples))
                        if all_clusters not in changes.keys():
                            changes[all_clusters] = []
                            changes[all_clusters].append(change_info)
                        else:
                            changes_all_clusters = changes[all_clusters][:]
                            for change in changes_all_clusters:
                                old,old_l,modification,new,new_l,n_new_samples,new_samples = change
                                if old == all_clusters and new == new_cluster_name:
                                    if modification == "new (merge)":
                                        changes[all_clusters].remove(change)
                                        changes[all_clusters].append(change_info)
                                else:
                                    changes[all_clusters].append(change_info)
                    elif len(exclusive_old) == 0 and len(exclusive_new) > 0: # cluster did not split and did merge
                        cluster_old_mx = mx_new[mx_new[mx_new.columns[0]].isin(list(cluster_samples))]
                        all_clusters = sorted(cluster_old_mx[old_partition].unique())
                        all_clusters = "_".join(all_clusters)
                        if all_clusters not in merged_clusters.keys():
                            max_cluster += 1
                            new_cluster_name = "cluster_" + str(max_cluster)
                            merged_clusters[all_clusters] = new_cluster_name
                        else:
                            new_cluster_name = merged_clusters[all_clusters]
                        conversion[cluster] = new_cluster_name
                        change_info = all_clusters,"multiple clusters","new (merge)",new_cluster_name,str(len(cluster_samples)),str(len(new_samples)),",".join(list(new_samples))
                        if all_clusters not in changes.keys():
                            changes[all_clusters] = []
                            changes[all_clusters].append(change_info)
                        else:
                            for change in changes[all_clusters]:
                                old,old_l,modification,new,new_l,n_new_samples,new_samples = change
                                if old == all_clusters and new != new_cluster_name:
                                    changes[all_clusters].append(change_info)
        elif len(new_samples) > 0: # cluster has new samples
            db_samples_lst = list(db_samples.keys())
            db_samples_lst.sort()
            for sample_list in db_samples_lst:
                old_samples = sample_list.split(",")
                old_cluster = str(db_samples[sample_list])
                cluster_intersection = set(cluster_samples) & set(old_samples)
                if len(cluster_intersection) > 0: # they share something
                    exclusive_old = set(old_samples) - set(cluster_samples)
                    exclusive_new = set(cluster_samples) - set(old_samples)
                    if len(exclusive_old) == 0 and len(exclusive_new) == len(new_samples): # cluster remains exactly the same and added new samples
                        if len(old_samples) == 1: # it was a singleton and now became a cluster
                            if "singleton_" in old_cluster:
                                max_cluster += 1
                                new_cluster_name = "cluster_" + str(max_cluster)
                                change_info = old_cluster,str(len(old_samples)),"new (increase)",new_cluster_name,str(len(cluster_samples)),str(len(new_samples)),",".join(list(new_samples))
                            else:
                                new_cluster_name = old_cluster
                                change_info = old_cluster,str(len(old_samples)),"kept (increase)",new_cluster_name,str(len(cluster_samples)),str(len(new_samples)),",".join(list(new_samples))
                            conversion[cluster] = new_cluster_name
                        else: # it was already a cluster that increased
                            conversion[cluster] = old_cluster
                            change_info = old_cluster,str(len(old_samples)),"kept (increase)",old_cluster,str(len(cluster_samples)),str(len(new_samples)),",".join(list(new_samples))
                        if old_cluster not in changes.keys():
                            changes[old_cluster] = []
                        changes[old_cluster].append(change_info)
                    elif len(exclusive_old) > 0 and len(exclusive_new) == len(new_samples): # cluster was split, increased and did not merge
                        if old_cluster not in split_clusters.keys():
                            split_clusters[old_cluster] = 0
                        if old_cluster not in changes.keys():
                            changes[old_cluster] = []
                        split_clusters[old_cluster] += 1
                        new_cluster_name = old_cluster + "." + str(split_clusters[old_cluster])
                        conversion[cluster] = new_cluster_name
                        change_info = old_cluster,str(len(old_samples)),"new (split_increase)",new_cluster_name,str(len(cluster_samples)),str(len(new_samples)),",".join(list(new_samples))
                        changes[old_cluster].append(change_info)
                    elif len(exclusive_old) > 0 and len(exclusive_new) > len(new_samples): # cluster was split, increase and did merge
                        cluster_old_mx = mx_new[mx_new[mx_new.columns[0]].isin(list(cluster_samples))]
                        all_clusters = sorted(cluster_old_mx[old_partition].unique())
                        all_clusters = "_".join(all_clusters)
                        if all_clusters not in merged_clusters.keys():
                            max_cluster += 1
                            new_cluster_name = "cluster_" + str(max_cluster)
                            merged_clusters[all_clusters] = new_cluster_name
                        else:
                            new_cluster_name = merged_clusters[all_clusters]
                        conversion[cluster] = new_cluster_name
                        change_info = all_clusters,"multiple clusters","new (split_merge_increase)",new_cluster_name,str(len(cluster_samples)),str(len(new_samples)),",".join(list(new_samples))
                        if all_clusters not in changes.keys():
                            changes[all_clusters] = []
                            changes[all_clusters].append(change_info)
                        else:
                            changes_all_clusters = changes[all_clusters][:]
                            for change in changes_all_clusters:
                                old,old_l,modification,new,new_l,n_new_samples,new_samples_name = change
                                if old == all_clusters and new == new_cluster_name:
                                    if modification == "new (merge_increase)":
                                        changes[all_clusters].remove(change)
                                        changes[all_clusters].append(change_info)
                                else:
                                    changes[all_clusters].append(change_info)
                    elif len(exclusive_old) == 0 and len(exclusive_new) > len(new_samples): # cluster did not split, increase and did merge
                        cluster_old_mx = mx_new[mx_new[mx_new.columns[0]].isin(list(cluster_samples))]
                        all_clusters = sorted(cluster_old_mx[old_partition].unique())
                        all_clusters = "_".join(all_clusters)
                        if all_clusters not in merged_clusters.keys():
                            max_cluster += 1
                            new_cluster_name = "cluster_" + str(max_cluster)
                            merged_clusters[all_clusters] = new_cluster_name
                        else:
                            new_cluster_name = merged_clusters[all_clusters]
                        conversion[cluster] = new_cluster_name
                        change_info = all_clusters,"multiple clusters","new (merge_increase)",new_cluster_name,str(len(cluster_samples)),str(len(new_samples)),",".join(list(new_samples))
                        if all_clusters not in changes.keys():
                            changes[all_clusters] = []
                            changes[all_clusters].append(change_info)
                        else:
                            for change in changes[all_clusters]:
                                old,old_l,modification,new,new_l,n_new_samples,new_samples_name = change
                                if old == all_clusters and new != new_cluster_name:
                                    changes[all_clusters].append(change_info)
        if cluster not in conversion.keys():
            if len(cluster_samples) == 1:
                max_singleton += 1
                new_cluster_name = "singleton_" + str(max_singleton)
                conversion[cluster] = new_cluster_name
            else:
                max_cluster += 1
                new_cluster_name = "cluster_" + str(max_cluster)
                conversion[cluster] = new_cluster_name

    return conversion, changes

def replace_cluster_names(conversion, mx):
    """ Replace the cluster names in the matrix
    input: dictionary with conversion names and data series
    output: new data series """
    
    mx.replace(conversion, inplace=True)
      
    return mx

def run_nomenclature(partitions, nomenclature, tag, day, log_name):
	""" use a previous nomenclature to replace cluster names in the partitions file 
	input: current partitions and nomenclature partitions in tsv
	output: partitions with new names as dataframe and tsv file """

	mx_original = pandas.read_table(partitions, dtype = str)
	mx_new = pandas.read_table(nomenclature, dtype = str)
	mx_new_flt, old_samples = filter_db(mx_original, mx_new)

	db_samples = get_new_db(mx_new)
    
	conversion = {} # conversion[partition][cluster in partitions] = cluster nomenclature
	changes_output = open(tag + "_nomenclature_changes.tsv", "w+")
	
	print("partition\told_cluster\told_cluster_length\tnomenclature_change\tcluster" + "_" + str(day) + "\tcluster_" + str(day) + "_length\tn_increase\tsamples_increase", file = changes_output)

	for partition in db_samples.keys():
		if partition in mx_original.columns[1:]: # exact match so only the column corresponding to the partition will be modified
			max_cluster, max_singleton = get_last_counters(mx_new[partition])
			conversion, changes = determine_new_cluster_names(old_samples, db_samples[partition], partition, mx_new_flt, partition, mx_original, max_cluster, max_singleton)
			mx_original[partition] = replace_cluster_names(conversion, mx_original[partition])
			if len(changes.keys()) > 0:
				out_partition = partition
				for cluster in sorted(changes.keys()):
					for old,old_l,modification,new,new_l,n_new_samples,new_samples in changes[cluster]:
						print(str(out_partition) + "\t" + str(old) + "\t" + str(old_l) + "\t" + str(modification) + "\t" + str(new) + "\t" + new_l + "\t" + n_new_samples + "\t" + new_samples, file = changes_output)
		else:
			if len(db_samples.keys()) == 1:
				print_log("\tThe nomenclature column " + str(partition) + " was not found in the partitions table. Cluster names will be modified for all possible partitions!", log_name)
				max_cluster, max_singleton = get_last_counters(mx_new[partition])
				for new_partition in mx_original.columns[1:]:
					conversion, changes = determine_new_cluster_names(old_samples, db_samples[partition], partition, mx_new_flt, new_partition, mx_original, max_cluster, max_singleton)
					mx_original[new_partition] = replace_cluster_names(conversion, mx_original[new_partition])
					if len(changes.keys()) > 0:
						out_partition = new_partition
						for cluster in sorted(changes.keys()):
							for old,old_l,modification,new,new_l,n_new_samples,new_samples in changes[cluster]:
								print(str(out_partition) + "\t" + str(old) + "\t" + str(old_l) + "\t" + str(modification) + "\t" + str(new) + "\t" + new_l + "\t", n_new_samples + "\t" + new_samples, file = changes_output)
			else:
				print_log("\tThe nomenclature column " + str(partition) + " was not found in the partitions table. Cluster names will not be modified!", log_name)
				changes = {}

	mx_original.to_csv(tag + "_partitions.tsv", index = False, header=True, sep ="\t")
	changes_output.close()

	return mx_original

def nomenclature_change2summary(tag):
	""" add column with nomenclature change to partitions_summary """

	partitions_mx = pandas.read_table(tag + "_partitions_summary.tsv", dtype = str)
	changes_mx = pandas.read_table(tag + "_nomenclature_changes.tsv", dtype = str)

	partitions_mx.set_index(["partition", "cluster"], drop = True, inplace = True)
	changes_mx.set_index(["partition", changes_mx.columns[4]], drop = True, inplace = True)
	
	new_mx = pandas.concat([changes_mx["nomenclature_change"].reindex(partitions_mx.index), changes_mx["n_increase"].reindex(partitions_mx.index), partitions_mx, changes_mx["samples_increase"].reindex(partitions_mx.index)], axis=1)

	new_mx.to_csv(tag + "_partitions_summary.tsv", index = True, header=True, sep = "\t")

def nomenclature_code(metadata_mx,metadata_col,partitions,tag,codes,day,log):
	""" add column with a hierarchical code to metadata_w_partitions """

	partitions_mx = pandas.read_table(partitions, dtype = str)
	partitions_mx.set_index(partitions_mx.columns[0], drop = False, inplace = True)
	metadata_mx.set_index(metadata_mx.columns[0], drop = False, inplace = True)
	original_columns = metadata_mx.columns.tolist()
	metadata_mx.columns = metadata_col

	col4nomenclature = codes
	list_codes = col4nomenclature[:]
	
	to_drop = []

	counter_partitions = 0
	counter_metadata = 0
	
	for code in list_codes:
		if code in partitions_mx.columns:
			if code not in metadata_col:
				metadata_mx = pandas.concat([metadata_mx, partitions_mx[code]], axis=1)
				to_drop.append(code)
				counter_partitions += 1
		elif code in metadata_col:
			counter_metadata += 1
			if counter_metadata > 1:
				col4nomenclature.remove(code)
				print_log("\tOnly one metadata code is allowed for the nomenclature code. So, " + code + " will be ignored!", log)
		else:
			print_log("\t" + str(code) + " was not found. So, it will be ignored!", log)
			col4nomenclature.remove(code)

	print_log("The nomenclature code has " + str(counter_partitions) + " partition levels and " + str(counter_metadata) + " metadata variable.", log)
	metadata_mx["nomenclature_code_" + str(day)] = metadata_mx[col4nomenclature].apply(lambda row: "-".join(row.values.astype(str)), axis=1)
	metadata_mx["nomenclature_code_" + str(day)].replace(" ", "_", regex=True, inplace=True)
	metadata_mx["nomenclature_code_" + str(day)].replace("cluster_", "C", regex=True, inplace=True)
	metadata_mx["nomenclature_code_" + str(day)].replace("singleton_", "S", regex=True, inplace=True)
	metadata_mx["nomenclature_code_" + str(day)].replace("-nan", "", regex=True, inplace=True)
	metadata_mx["nomenclature_code_" + str(day)].replace("nan", "", regex=True, inplace=True)

	for c in set(to_drop):
		metadata_mx.drop(c, axis=1, inplace=True)  
	original_columns.append("nomenclature_code_" + str(day))
	metadata_mx.columns = original_columns 
	metadata_mx.to_csv(tag + "_metadata_w_partitions.tsv", index = False, header=True, sep ="\t")

def solve_nomenclature_code_in_args(args,day):
	""" solve args inputs when the user requests the nomenclature_code """
	
	# metadata2report
	if args.metadata2report != "none":
		if "nomenclature_code" in args.metadata2report:
			metadata2report = str(args.metadata2report).split(",")
			for i in range(len(metadata2report)):
				if metadata2report[i] == "nomenclature_code":
					metadata2report[i] = "nomenclature_code_" + str(day)
			args.metadata2report = ",".join(metadata2report)
	else:
		args.metadata2report = "nomenclature_code_" + str(day)
	
	# zoom-in
	if args.zoom != "no" and "nomenclature_code" in args.zoom:
		zoom = str(args.zoom).split(",")
		for i in range(len(zoom)):
			if zoom[i] == "nomenclature_code":
				zoom[i] = "nomenclature_code_" + str(day)
		args.zoom = ",".join(zoom)
	
	# frequency matrix
	if args.frequency_matrix != "no":
		final_requests = []
		for combination in args.frequency_matrix.split(";"):
			v1,v2 = combination.split(",")
			if v1 == "nomenclature_code":
				v1 = "nomenclature_code_" + str(day)
			elif v2 == "nomenclature_code":
				v2 = "nomenclature_code_" + str(day)
			request = str(v1) + "," + str(v2)
			final_requests.append(request)
		args.frequency_matrix = ";".join(final_requests)

	# count matrix
	if args.count_matrix != "no":
		final_requests = []
		for combination in args.count_matrix.split(";"):
			v1,v2 = combination.split(",")
			if v1 == "nomenclature_code":
				v1 = "nomenclature_code_" + str(day)
			elif v2 == "nomenclature_code":
				v2 = "nomenclature_code_" + str(day)
			request = str(v1) + "," + str(v2)
			final_requests.append(request)
		args.count_matrix = ";".join(final_requests)

	return args

def rename_clusters_subsets(analysis,partitions,tag):
	""" rename clusters resulting from zoom-in or subtree analysis """

	mx = pandas.read_table(partitions)

	new_name = {}
	if analysis == "zoom-in":
		cluster_zoom_in = "_".join(tag.split("_")[-2:])
		new_name["cluster_"] = cluster_zoom_in + "_group_"
		new_name["singleton_"] = cluster_zoom_in + "_singleton_"
	elif analysis == "subtree":
		code = tag.split("_")[-1]
		new_name["cluster_"] = code + "_group_"
		new_name["singleton_"] = code + "_singleton_"
	
	mx.replace(new_name, regex=True, inplace=True)
	mx.to_csv(partitions, index = False, header=True, sep = "\t")


# running the pipeline	----------

def main():
    
	# argument options	----------
    
	parser = argparse.ArgumentParser(prog="reportree.py", formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent("""\
									###############################################################################             
									#                                                                             #
									#                                reportree.py                                 #
									#                                                                             #
									############################################################################### 
									                            
									ReporTree aims to answer the need of:
									
									-obtaining genetic clusters at any distance thresholds of a rooted tree (e.g. 
									SNP-scaled tree), SNP or cg/wgMLST allele matrix, VCF files, sequence 
									alignment, or distance matrix

									-obtaining summary reports with the statistics/trends (e.g. timespan, location,
									cluster/group size and composition, age distribution etc.) for the derived 
									genetic clusters or for any other provided grouping variable (such as, clade, 
									lineage, ST, vaccination status, etc.)
									
									-obtaining count/frequency matrices for the derived genetic clusters or for any 
									other provided grouping variable
									
									-identifying regions of cluster stability (i.e. partition threshold ranges in 
									which cluster composition is similar)
									
									
									WARNING 'white spaces':
									- they should be avoided in metadata and partition tables column 
									names (if you include them, please use an underscore in the command line 
									instead of the white space)
									- White cells in the metadata body, they will be replaced by 'EMPTY'
									- they should also be avoided in sample names because we will only consider the 
									first part of the name.
									
									Note 2: To take the most profit of this script we recommend that you include 
									the column 'date' in the metadata. This column must follow the format YYYY-MM-DD
									If you only provide YYYY, it will assume YYYY-01-01!!
									
									Note 3: If a 'date' column is provided in the metadata, this script will 
									copy it to the "date_original" column (which remains intact) and create a new
									"date" column with the date formated as YYYY-MM-DD. A new "year" column will
									also be created, corresponding to the exact year of the sample (not ISO). 
									Moreover, it will determine and provide in the new metadata table the columns 
									iso_year, iso_week_nr and iso_week for each sequence (e.g. iso_year = 2021, 
									iso_week_nr = 52, iso_week = 2021W52)!!
									
									Note 4: While for nominal or categorical variables this script can provide in 
									the summary report the number of observations or the frequency of each 
									observation, for the 'date' column this script can provide:
										- first_seq_date
										- last_seq_date
										- timespan_days
									
									Note 5: Currently, for non-SNP-distance rooted trees, users must specify a 
									minimum unit to cut the tree (currently, the default is 1, which is equivalent 
									to 1 SNP in a SNP-scaled rooted tree).
									
									Note 6: By default, we assume '0' as the missing data code except in sequence 
									alignments in which we assume 'N' as missing data code. If your code differs from 
									this, please use the '--missing-code' argument.									
									
									TIP!! If you do not know which columns you can indicate for the argument 
									'--columns_summary_report', use the '--list' argument!! We strongly advise you to
									run this argument in the first time you are running ReporTree :-)
									
									
									Check the github page for information about citation.
									                  
									-------------------------------------------------------------------------------"""))
	
	## general parameters
	
	versioning = parser.add_argument_group("Version", "ReporTree version")
	versioning.add_argument("-v", "--version", dest="version", action="store_true", help="Print version and exit")
	
	group0 = parser.add_argument_group("ReporTree", "ReporTree input/output specifications")
	group0.add_argument("-a", "--allele-profile", dest="allele_profile", default="", required=False, type=str, help="Input allele/SNP profile matrix (tsv format)")
	group0.add_argument("-align", "--alignment", dest="alignment", default="", required=False, type=str, help="Input multiple sequence alignment (fasta format)")
	group0.add_argument("-d_mx", "--distance_matrix", dest="distance_matrix", default="", required=False, type=str, help="Input pairwise distance matrix (tsv format)")
	group0.add_argument("-t", "--tree", dest="tree", default="", required=False, type=str, help="Input tree (newick format)")
	group0.add_argument("-p", "--partitions", dest="partitions", required=False, default="", type=str, help="Partitions file (tsv format) - 'partition' represents the threshold at \
						which clustering information was obtained")
	group0.add_argument("-m", "--metadata", dest="metadata", required=False, type=str, default = "none", help="[MANDATORY] Metadata file (tsv format). To take the most profit of ReporTree \
						functionalities, you must provide this file.")
	group0.add_argument("-vcf", "--vcf", dest="vcf", default="", required=False, type=str, help="Single-column list of VCF files (txt format). This file must comprise the full PATH to \
						each vcf file.")
	group0.add_argument("-var", "--variants", dest="variants", default="", required=False, type=str, help="Input table (tsv format) with sample name in the first column and a \
						comma-separated list of variants in the second column with the following regular expression: '\w(\d+)\w' ")
	group0.add_argument("-out", "--output", dest="output", required=False, default="ReporTree", type=str, help="[OPTIONAL] Tag for output file name (default = ReporTree)")
	group0.add_argument("--list", dest="list_col_summary", required=False, action="store_true", help="[OPTIONAL] If after your command line you specify this option, ReporTree will list all the \
						possible columns that you can use as input in '--columns_summary_report'. To obtain information about the partition name for other arguments('--frequency-matrix' and/or \
		     			'--count-matrix'), please also indicate the type of analysis. NOTE!! The objective of this argument is to help you with the input of some other arguments. So, it will not \
		     			run reportree.py main functions!!")
	
	
	## analysis details
	
	group1 = parser.add_argument_group("Analysis details", "Analysis details")
	group1.add_argument("--analysis", dest="analysis", required=False, default="", type=str, help="Type of clustering analysis (options: grapetree, HC, treecluster). If you provide a tree, \
						genetic clusters will always be obtained with treecluster. If you provide a distance matrix, genetic clusters will always be obtained with HC. If you provide any other \
						input, it is MANDATORY to specify this argument.")
	group1.add_argument("--subset", dest="subset", required=False, action="store_true", help="[OPTIONAL] Obtain genetic clusters using only the samples that correspond to the filters specified in the \
						'--filter' argument.")
	group1.add_argument("-d", "--dist", dest="dist", required=False, default=1.0, type=float, help="[OPTIONAL] Distance unit by which partition thresholds will be multiplied (example: if -d 10 and \
						-thr 5,8,10-30, the minimum spanning tree will be cut at 50,80,100,110,120,...,300. If -d 10 and --method-threshold avg_clade-2, the avg_clade threshold will be set \
						at 20). This argument is particularly useful for non-SNP-scaled trees. Currently, the default is 1, which is equivalent to 1 allele distance or 1 SNP distance. [1.0]")
						
			
	## cleaning missing data
	
	group2 = parser.add_argument_group("Cleaning missing data", "Remove loci/positions and samples based on missing content")
	group2.add_argument("--missing-code", dest="missing_code", required=False, type=str, default = "0", help="[OPTIONAL] Code representing missing data. If different from '0' or 'N', please \
		     			try to avoid a IUPAC character (even in lower-case). [default: 'N' when '-align' provided and '0' for other inputs]")
	group2.add_argument("--site-inclusion", dest="N_content", required=False, default = 0.0, help="[OPTIONAL: Useful to remove informative sites/loci with excess of missing data] Minimum \
						proportion of samples per site without missing data (e.g. '--site-inclusion 1.0' will only keep loci/positions without missing data, i.e. a core alignment/profile; \
						'--site-inclusion 0.0' will keep all loci/positions) NOTE: This argument works on profile/alignment loci/positions (i.e. columns)! [default: 0.0].")
	group2.add_argument("--loci-called", dest="loci_called", required=False, default = 0.0, help="[OPTIONAL - only works for matrices] Minimum proportion of loci/positions called for SNP/allele \
		     			matrices (e.g. '--loci-called 0.95' will only keep in the profile matrix samples with > 95%% of alleles/positions, i.e. <= 5%% missing data). Applied after '--site-inclusion' \
						argument! [default: 0.0]")
	group2.add_argument("--sample-ATCG-content", dest="ATCG_content", required=False, default = 0.0, help="[OPTIONAL - only works for alignment] Minimum proportion (0 to 1) of ATCG in informative \
		     			sites of the alignment per sample (e.g. '--sample-ATCG-content 1.0' will only keep samples without N's or any non-ATCG code in informative sites) [default: 0 - keep all samples]")
	
	
	# alignment processing
	
	group3 = parser.add_argument_group("Alignment processing", "Alignment processing")
	group3.add_argument("--remove-reference", dest="remove_ref", required=False, action="store_true", help="Set only if you want to remove the reference sequence of the alignment (reference \
						name must be provided with the argument '--reference').")					
	group3.add_argument("--use-reference-coords", dest="use_ref", required=False, action="store_true", help="Set only if you want that column names in the final alignment matrix represent the \
						reference coordinates (reference name must be provided with the argument '--reference') Note: Depending on the alignment size, this argument can make alignment processing \
						very slow!")
	group3.add_argument("-r", "--reference", dest="reference", required=False, type=str, default = "none", help="[OPTIONAL] Name of reference sequence. Required if '--remove-reference' and/or \
						'--use-reference-coords' specified.")
	group3.add_argument("--get-position-correspondence", dest="pos_corr", required=False, default="none", help="[OPTIONAL] Request a .tsv with position correspondence between any sequences \
						of your alignment. These should be indicated separated by a comma (e.g. seqA,seqB). To get the position coordinates of all sequences just write 'all'.")
	group3.add_argument("--position-list", dest="pos_list", required=False, default="none", help="[OPTIONAL] .tsv file with the positions of interest to be reported when \
						'--get-position-correspondence' is requested. Each column should correspond to the positions of a sequence and the sequence name should be indicated in the header. If this \
						file is not provided, all positions of the alignment will be reported.")
	
	
	## partitioning grapetree
	
	group4 = parser.add_argument_group("Partitioning with GrapeTree", "Specifications to get and cut minimum spanning trees")
	group4.add_argument("--method", dest="grapetree_method", default="MSTreeV2", help="\"MSTreeV2\" [DEFAULT]\n Alternative:\"MSTree (goeBURST)\"\n")
	group4.add_argument("--missing", dest="handler", default=0, type=int, help="ONLY FOR MSTree. \n0: [DEFAULT] ignore missing data in pairwise comparison. \n1: remove column \
						with missing data. \n2: treat missing data as an allele. \n3: use absolute number of allelic differences.")
	group4.add_argument("--n_proc",  dest="number_of_processes", type=int, default=5, help="Number of CPU processes in parallel use. [5]")
	group4.add_argument("-thr", "--threshold", dest="threshold", default = "max", help="[OPTIONAL] Partition threshold for clustering definition (integer). Different thresholds can be comma-separated (e.g. 5,8,16). \
						Ranges can be specified with a hyphen separating minimum and maximum (e.g. 5,8,10-20). If this option is not set, the script will perform clustering for all the values in the range 0 to max. \
		     			If you prefer to exclusively use the '-pct_thr' argument, please set '-thr none'. Note: Threshold values are inclusive, i.e. '-thr 7' will consider samples with <= 7 differences as belonging \
		     			to the same cluster!")
	group4.add_argument("-pct_thr", "--pct_threshold", dest="pct_threshold", default = "none", help="[OPTIONAL] Similar to 'thr' but values are indicated as the proportion of differences to the \
						final allelic schema size or number of informative positions, e.g. '-pct_thr 0.005' corresponds to a threshold of 5 allelic/SNP differences in a matrix with 1000 \
						loci/sites under analysis). Different values can be comma-separated (e.g. 0.005,0.01,0.1). Ranges CANNOT be specified. This option is particularly useful for dynamic \
						wgMLST analysis for which the size of the schema under analysis is contigent on dataset diversity. Note: This argument can be specified even if you used the '-thr' \
						argument.")
	group4.add_argument("--matrix-4-grapetree", dest="matrix4grapetree", required=False, action="store_true", help="Output an additional allele profile matrix with the header ready for GrapeTree \
						visualization. Set only if you WANT the file!")
	group4.add_argument("--wgMLST", dest="wgmlst", default=False, action="store_true", help="[EXPERIMENTAL] a better support of wgMLST schemes (check GrapeTree github for details).")
	
	
	## partitioning hierarchical clustering
	
	group5 = parser.add_argument_group("Partitioning with HC", "Specifications to genetic clusters with hierarchical clustering")
	group5.add_argument("--HC-threshold", dest="HCmethod_threshold", required=False, default="single", help="[OPTIONAL] List of HC methods and thresholds to include in the analysis (comma-separated). To \
						get clustering at all possible thresholds for a given method, write the method name (e.g. single). To get clustering at a specific threshold, indicate the threshold with \
						a hyphen (e.g. single-10). To get clustering at a specific range, indicate the range with a hyphen separating minimum and maximum (e.g. single-2-10). If you prefer to exclusively use the \
						'--pct-HC-threshold' argument, please set '--HC-threshold none'. Note: Threshold values are inclusive, i.e. '--HC-threshold single-7' will consider samples with <= 7 \
						differences as belonging to the same cluster! Default: single (Possible methods: single, complete, average, weighted, centroid, median, ward)")
	group5.add_argument("--pct-HC-threshold", dest="pct_HCmethod_threshold", required=False, default="none", help="[OPTIONAL] Similar to '--HC-threshold' but the partition threshold for cluster definition \
						is set as the proportion of differences to the final allelic schema size or number of informative positions, e.g. '--pct-HC-threshold single-0.005' corresponds to a \
						threshold of 5 allelic/SNP differences in a matrix with 1000 loci/sites under analysis. Ranges CANNOT be specified.")
											
						
	## partitioning treecluster
	
	group6 = parser.add_argument_group("Partitioning with TreeCluster", "Specifications to cut the tree with TreeCluster")
	group6.add_argument("--method-threshold", dest="method_threshold", required=False, default="root_dist,avg_clade-1", 
						help="List of TreeCluster methods and thresholds to include in the analysis (comma-separated). To get clustering at all possible thresholds for a given method, write \
						the method name (e.g. root_dist). To get clustering at a specific threshold, indicate the threshold with a hyphen (e.g. root_dist-10). To get clustering at a specific \
						range, indicate the range with a hyphen (e.g. root_dist-2-10). Default: root_dist,avg_clade-1 (List of possible methods: avg_clade, leaf_dist_max, leaf_dist_min, length, \
						length_clade, max, max_clade, root_dist, single_linkage, single_linkage_cut, single_linkage_union) Warning!! So far, ReporTree was only tested with avg_clade and \
						root_dist!")
	group6.add_argument("--support", dest="support", required=False, default="-inf", help="[OPTIONAL: see TreeCluster github for details] Branch support threshold") 
	group6.add_argument("--root-dist-by-node", dest="root_dist_by_node", required=False, action="store_false", help="[OPTIONAL] Set only if you WANT to cut the tree with root_dist method at each tree \
						node distance to the root (similar to root_dist at all levels but just for informative distances)")
	group6.add_argument("-root", "--root", dest="root", required=False, default="no", help="Set root of the input tree. Specify the leaf name to use as output. Alternatively, write 'midpoint', \
						if you want to apply midpoint rooting method.")
	
	
	## cluster nomenclature
	group7 = parser.add_argument_group("ReporTree cluster nomenclature", "Cluster nomenclature instructions")
	group7.add_argument("--nomenclature-file", dest="nomenclature", required=False, type=str, default="", help="[OPTIONAL] Intended usage: provide a .tsv file with the nomenclature information to be \
		     			used (normaly the '*_partitions.tsv' file of a previous ReporTree run) with the following structure: First column should have the sample names; Subsequent columns have the partitions \
		     			at any/all levels. Columns matching the column names of the partitions requested in the current run (e.g. MST-7x1.0) will be used to name the current clusters at the respective \
		     			levels. For more details on how to use this nomenclature system please visit our github! Alternative usage (at your own risk): you can provide a .tsv file with just a single column with \
		     			a grouping variable that does not match any partition requested in the current run (e.g. ST). In this case, all clusters will be named according to this column (e.g. ST1.1,ST1.2,ST2.1...).")
	group7.add_argument("--nomenclature-code-levels", dest="code_levels", required=False, type=str, default="", help="[OPTIONAL and only available for --analysis grapetree or HC] This argument allows getting a \
		     			nomenclature code combining cluster information at different hierarchical levels (e.g. if '150,30,7' is provided when the analysis is GrapeTree, a code combining the cluster names at these \
		     			levels will be generated: C3-C2-C1). The order of levels indicated in the command-line will be kept in the code (in the example, C3 indicates that the sample belongs to cluster_3 at MST-150). \
		     			Of note, if a '--nomenclature-file' is provided in subsequent ReporTree runs using the same method and thresholds, the nomenclature code will be kept. You can also add one metadata variable \
		     			(e.g. Country) to get an extra layer to the code (e.g. C3-C2-C1-Portugal). Partition thresholds can be indicated in this argument following the same rules as the arguments '-thr' and '-pct_thr' \
		     			for GrapeTree or '--HC-threshold' and '--pct-HC-threshold' for HC.")


	## reportree
	
	group8 = parser.add_argument_group("ReporTree metadata report", "Specific parameters to report clustering/grouping information associated to metadata")
	group8.add_argument("--columns_summary_report", dest="columns_summary_report", required=False, default="n_sequence,lineage,n_country,country,n_region,first_seq_date,last_seq_date,timespan_days", type=str, help="Columns \
						(i.e. variables of metadata) to get statistics for the derived genetic clusters or for other grouping variables defined in --metadata2report (comma-separated). If the \
						name of the column is provided, the different observations and the respective percentage are reported. If 'n_column' is specified, the number of the different \
						observations is reported. For example, if 'n_country' and 'country'  are specified, the summary will report the number of countries and their distribution (percentage) \
						per cluster/group. Exception: if a 'date' column is in the metadata, it can also report first_seq_date, last_seq_date, timespan_days. Check '--list' argument for some help. \
						Default = n_sequence,lineage,n_country,country,n_region,first_seq_date,last_seq_date,timespan_days [the order of the list will be the order of the columns in the report]")
	group8.add_argument("--partitions2report", dest="partitions2report", required=False, default="all", type=str, help="Thresholds for which clustering information will be included in a joint report \
						(comma-separated). Other alternatives: 'all' == all partitions; 'stability_regions' == first partition of each stability region as determined by \
						comparing_partitions_v2.py. Note: 'stability_regions' can only be inferred when partitioning TreeCluster or GrapeTree or HC is run for all possible thresholds or when a \
						similar partitions table is provided (i.e. sequential partitions obtained with the same clustering method) [all]. Partition thresholds can be indicated in this argument \
		     			following the same rules as the arguments '-thr' and '-pct_thr' for GrapeTree or '--HC-threshold' and '--pct-HC-threshold' for HC or '--method-threshold' for TreeCluster.")
	group8.add_argument("--metadata2report", dest="metadata2report", required=False, default="none", help="Columns of the metadata table for which a separated summary report must be provided \
						(comma-separated)")
	group8.add_argument("-f", "--filter", dest="filter_column", required=False, default="", help="[OPTIONAL] Filter for metadata columns to select the samples to analyze. This must be specified \
						within quotation marks in the following format 'column< >operation< >condition' (e.g. 'country == Portugal'). When more than one condition is specified for a given column, \
						they must be separated with commas (e.g 'country == Portugal,Spain,France'). When filters include more than one column, they must be separated with semicolon (e.g. \
						'country == Portugal,Spain,France;date > 2018-01-01;date < 2022-01-01'). White spaces are important in this argument, so, do not leave spaces before and after \
						commas/semicolons.")
	group8.add_argument("--sample_of_interest", dest="sample_of_interest", required=False, default="all", help="[OPTIONAL] List of samples of interest for which summary reports will be created. This \
		     			list can be a comma-separated list in the command line, or a comma-separated list in a file, or a list in the first column of a tsv file. No headers should be provided in the \
		     			input files. If nothing is provided, only the summary reports comprising all samples will be generated.")
	group8.add_argument("--zoom-cluster-of-interest", dest="zoom", required=False, default="no", help="[OPTIONAL and only available for --analysis grapetree or HC] Repeat the analysis using only the \
		     			samples that belong to each cluster of the samples of interest at a given distance threshold. This argument takes as input a comma-separated list of partitions for which you \
		     			want the zoom-in. Partition thresholds can be indicated in this argument following the same rules as the arguments '-thr' and '-pct_thr' for GrapeTree or '--HC-threshold' and \
		     			'--pct-HC-threshold' for HC. This argument requires that a metadata table was provided with '-m'. Default: no zoom-in.")
	group8.add_argument("--subtree-of-interest", dest="subtree", required=False, default="no", help="[OPTIONAL and only available for --analysis grapetree or HC] Repeat the analysis using the n \
		     			closest samples of each sample of interest. This argument takes as input a comma-separated list of n's, corresponding to the number of closest samples you want to include for \
		     			the samples of interest. This argument requires that a metadata table was provided with '-m'. Default: no subtree.")
	group8.add_argument("--unzip", dest="unzip", required=False, action="store_true", help="[OPTIONAL and only available for --analysis grapetree or HC] Provide the outputs of '--zoom-cluster-of-interest' and \
						'--subtree-of-interest' in unzipped format.")
	group8.add_argument("--frequency-matrix", dest="frequency_matrix", required=False, default="no", help="[OPTIONAL] Metadata column names for which a frequency matrix will be generated. This \
						must be specified within quotation marks in the following format 'variable1,variable2'. Variable1 is the variable for which frequencies will be calculated (e.g. for \
						'lineage,iso_week' the matrix reports the percentage of samples that correspond to each lineage per iso_week). If you want more than one matrix you can separate the \
						different requests with semicolon (e.g. 'lineage,iso_week;country,lineage'). If you want a higher detail in your variable2 and decompose it into two columns you use a \
						colon (e.g. lineage,country:iso_week will report the percentage of samples that correspond to each lineage per iso_week in each country)")
	group8.add_argument("--count-matrix", dest="count_matrix", required=False, default="no", help="[OPTIONAL] Same as '--frequency-matrix' but outputs counts and not frequencies")
	group8.add_argument("--mx-transpose", dest="mx_transpose", required=False, action="store_true", help="[OPTIONAL] Set ONLY if you want that the variable1 specified in '--frequency-matrix' \
						or in '--count-matrix' corresponds to the matrix first column.")
	group8.add_argument("--pivot", dest="pivot", required=False, action="store_true", help="[OPTIONAL] Set ONLY if you want an additional table for each count/frequency matrix in pivot format.")
						
						
	## comparing partitions
	
	group9 = parser.add_argument_group("Stability regions", "Congruence analysis of cluster composition at all possible partitions to determine regions of cluster stability (automatically run if you set --partitions2report 'stability_regions'). WARNING! This option is planned to handle sequential partitions obtained with the same clustering method, such as a partitions table derived from cg/wgMLST data (from 1 to max allele threshold). Use it at your own risk, if you provide your own partitions table.")
	group9.add_argument("-AdjW", "--AdjustedWallace", dest="AdjustedWallace", action= "store", default=0.99, help="Threshold of Adjusted Wallace score to consider an observation for method \
						stability analysis [0.99]")
	group9.add_argument("-n", "--n_obs", dest="n_obs", action="store", default=5, help="Minimum number of sequencial observations that pass the Adjusted Wallace score to be considered a \
						'stability region' (i.e. a threshold range in which cluster composition is similar) [5]")
	group9.add_argument("-o", "--order", dest="order", action= "store", default=0, required=False, help="[Set only if you provide your own partitions table] Partitions order in the partitions \
						table (0: min -> max; 1: max -> min) [0]")
	group9.add_argument("--keep-redundants", dest="keep_redundants", action= "store_true", help="Set ONLY if you want to keep all samples of each cluster of the most discriminatory partition\
						 (by default redundant samples are removed to avoid the influence of cluster size)")
	
						 
	args = parser.parse_args()


	# check if version	----------
	
	if args.version:
		print("version:", version, "\nlast_updated:", last_updated)
		sys.exit()
		
	
	# check if the user wants the list of columns	----------
	
	cmds = []
	if args.list_col_summary:
		if  args.metadata != "none":
			columns_metadata, columns_methods = col_list(args)
			print("\nPossible column names to provide as input to '--columns_summary_report' argument:")
			print("\n".join(columns_metadata))
		else:
			print("\nAs no metadata is provided, we do not know what you can indicate in the '--columns_summary_report' argument.")
		
		if args.analysis != "":
			print("\nReporTree will output some partition columns that can be included in your summary reports. For each method that you requested, the column name that you can indicate in '--columns_summary_report', '--frequency-matrix' and/or '--count-matrix' arguments must follow this structure:")
			for method in columns_methods:
				print(method + "-<threshold>x<dist>          ---------->          e.g. " + method + "-30x" + str(args.dist))
			print("\nFor the '--partitions2report' argument, if you did not provide a partitions table, you can indicate the partitions in the same way you indicated in the partition thresholds for the analysis.")
		print("\nThere are also some special codes that you can specify. For example:")
		print("-> if you write stability_regions in the '--partitions2report' argument, ReporTree will run comparing_partitions_v2.py and use the first partition of each stable region in this argument.")
		print("-> if you write nomenclature_code in the '--metadata2report', or the '--zoom-cluster-of-interest', or the '--frequency-matrix', or the '--count-matrix', ReporTree will infer what is the name of the column with the nomenclature_code of this run and provide the reports/metrics you requestes.")
		print("\n")
		sys.exit()
		
	
	# starting logs	----------
    
	log_name = args.output + ".log"
	log = open(log_name, "w+")

	print_log("\n******************** running reportree.py ********************\n", log)
	print_log("version " + str(version) + " last updated on " + str(last_updated) + "\n", log)
	print_log(" ".join(sys.argv), log)
	
	start = datetime.datetime.now()
	print_log("start: " + str(start), log)

	input_align = False
	day = str(datetime.datetime.today().strftime("%Y-%m-%d"))

	if args.metadata == "none":
		print_log("\nWARNING!! You did not provide a metadata file... metadata_report.py will not be run! Provide such a file if you want to take the most profit of ReporTree :-)", log)
		metadata_col = []
	else:
		metadata_mx = pandas.read_table(args.metadata)
		metadata_col = [col.replace(" ", "_") for col in metadata_mx.columns]
	
	# adapting inputs to report nomenclature_code
	if args.code_levels != "":
		args = solve_nomenclature_code_in_args(args,day)
	
    # reportree workflow	----------
    
    ## partitions table provided	--------------------
    
	if args.partitions != "": # partitions table provided
		if args.tree != "": # tree was provided -> need to get partitions again?
			print_log("\nTree and partitions files specified... I am confused :-(\n", log)
			sys.exit(1)
		
		elif args.allele_profile != "": # alleles were provided -> need to get partitions again?
			print_log("\nAllele profiles and partitions files specified... I am confused :-(\n", log)
			sys.exit(1)
		
		elif args.alignment != "": # alignment was provided -> need to get partitions again?
			print_log("\nSequence alignment and partitions files specified... I am confused :-(\n", log)
			sys.exit(1)
		
		elif args.vcf != "": # VCF was provided -> need to get partitions again?
			print_log("\nVCF and partitions files specified... I am confused :-(\n", log)
			sys.exit(1)
		
		elif args.variants != "": # variants was provided -> need to get partitions again?
			print_log("\nVariants list and partitions files specified... I am confused :-(\n", log)
			sys.exit(1)
		
		elif args.distance_matrix != "": # distance mx was provided -> need to get partitions again?
			print_log("\nDistance matrix and partitions files specified... I am confused :-(\n", log)
			sys.exit(1)
		
		else: # can continue using partitions
			print_log("\nPartitions file provided...\n", log)
			partition_status = "old"

			if args.nomenclature != "":
				print_log("\nNomenclature file provided. Will work on cluster names...", log)
				partitions = run_nomenclature(args.partitions, args.nomenclature, args.output, day, log)
				args.partitions = args.output + "_partitions.tsv"
				partition_status = "new"
				cmds.append("nomenclature")
			
			# running comparing partitions
			if "stability_regions" in args.partitions2report and "all" in args.partitions2report:
				print_log("\t'stability_regions' and 'all' options cannot be simultaneousl-y specified in --partitions2report... I am confused :-(", log)
				sys.exit(1)
			elif "stability_regions" in args.partitions2report and "all" not in args.partitions2report:
				print_log("\t'stability_regions' option specified. Will also run comparing_partitions_v2.py to determine the partitions that must be included in the report...", log)
				log.close()
				returned_value = run_stability(args, log_name, partition_status)
				if str(returned_value) != "0":
					sys.exit("\n\nReporTree exited before expected while running comparing_partitions_v2.py :-( Please check your input files and your command line. If you are in trouble and cannot figure out what is going on, contact us!!")
				log = open(log_name, "a+")
				
			partitions2report_final, partitions2report_lst = infer_partitions_from_list("other", args.partitions, metadata_col, args.partitions2report, args.output, args.dist, log)
			
			if len(partitions2report_final) == 0:
				print_log("\tThe analysis of the partitions to report returned an empty list. All partitions will be included in the report...", log)
				partitions2report_final = "all"
			
			# adding nomenclature_code to metadata

			if args.code_levels != "" and args.metadata != "none":
				print_log("\nGenerating nomenclature code...", log)
				code_levels, code_levels_lst = infer_partitions_from_list("other", args.partitions, metadata_col, args.code_levels, args.output, args.dist, log)
				nomenclature_code(metadata_mx,metadata_col,args.partitions,args.output,code_levels_lst,day,log)
				args.metadata = args.output + "_metadata_w_partitions.tsv"

			# getting metadata report
			if args.metadata != "none":
				log.close()
				returned_value = run_metadata_report(args, partitions2report_final, partition_status)
				if str(returned_value) != "0":
					sys.exit("\n\nReporTree exited before expected while running metadata_report.py :-( Please check your input files and your command line. If you are in trouble and cannot figure out what is going on, contact us!!")
				cmds.append("metadata")
				log = open(log_name, "a+")
			
			# adding nomenclature change to partitions summary
			if "nomenclature" in cmds and "metadata" in cmds:
				nomenclature_change2summary(args.output)
				
	
	## tree provided	--------------------
	
	elif args.tree != "": # tree was provided
		if args.allele_profile != "": # allelic profiles provided -> grapetree, HC or treecluster
			print_log("\nTree and allele profiles files specified... I am confused :-(\n", log)
			sys.exit(1)
		
		elif args.alignment != "": # alignment was provided -> grapetree, HC or treecluster
			print_log("\nTree and sequence alignment files specified... I am confused :-(\n", log)
			sys.exit(1)
		
		elif args.vcf != "": # VCF was provided -> grapetree, HC or treecluster
			print_log("\nTree and VCF files specified... I am confused :-(\n", log)
			sys.exit(1)
		
		elif args.variants != "": # variants was provided -> grapetree, HC or treecluster
			print_log("\nTree and variants list files specified... I am confused :-(\n", log)
			sys.exit(1)
		
		elif args.distance_matrix != "": # distance mx was provided -> grapetree, HC or treecluster
			print_log("\nTree and distance matrix specified... I am confused :-(\n", log)
			sys.exit(1)
		
		else: # can continue using tree and run treecluster and metadata report
			print_log("\nTree file provided -> will run partitioning_treecluster.py:\n", log)
			
			# running partitioning treecluster

			log.close()
			returned_value = run_treecluster(args)
			if str(returned_value) != "0":
				sys.exit("\n\nReporTree exited before expected while running partitioning_treecluster.py :-( Please check your input files and your command line. If you are in trouble and cannot figure out what is going on, contact us!!")
			log = open(log_name, "a+")

			partition_status = "new"
			if args.nomenclature != "":
				print_log("\nNomenclature file provided. Will work on cluster names...", log)
				partitions = run_nomenclature(args.output + "_partitions.tsv", args.nomenclature, args.output, day, log)
				partition_status = "new"
				cmds.append("nomenclature")
		
		# running comparing partitions
		if "stability_regions" in args.partitions2report and "all" in args.partitions2report:
			print_log("\t'stability_regions' and 'all' options cannot be simultaneously specified in --partitions2report... I am confused :-(", log)
			sys.exit(1)
		elif "stability_regions" in args.partitions2report and "all" not in args.partitions2report:
			if len(all_partitions_available(args.method_threshold, args.root_dist_by_node)) >= 1 and args.dist == 1.0: # at least all partitions for 1 method
				print_log("\t'stability_regions' option specified. Will also run comparing_partitions_v2.py to determine the partitions that must be included in the report...", log)
				for method in all_partitions_available(args.method_threshold, args.root_dist_by_node):
					filter_partitions_table(method, args.output + "_partitions.tsv", args.output)
					log.close()
					real_partitions = args.partitions
					args.partitions = args.output + "_tmp.tsv"
					real_output = args.output
					args.output = args.output + "_" + method
					returned_value = run_stability(args, log_name, partition_status)
					if str(returned_value) != "0":
						sys.exit("\n\nReporTree exited before expected while running comparing_partitions_v2.py :-( Please check your input files and your command line. If you are in trouble and cannot figure out what is going on, contact us!!")
					args.partitions = real_partitions
					args.output = real_output
					os.system("rm " + args.output + "_tmp.tsv")
					log = open(log_name, "a+")
			else:
				if len(all_partitions_available(args.method_threshold, args.root_dist_by_node)) >= 1 and args.dist != 1.0:
					print_log("\t'stability_regions' option specified but minimum distance was different from 1.", log)
					sys.exit(1)
				else:
					print_log("\t'stability_regions' option specified but no method was run for all possible partitions.", log)
					sys.exit(1)
			
		partitions2report_final, partitions2report_lst = infer_partitions_from_list("treecluster", args.output + "_partitions.tsv", metadata_col, args.partitions2report, args.output, args.dist, log)
		
		if len(partitions2report_final) == 0:
			print_log("\tThe analysis of the partitions to report returned an empty list. All partitions will be included in the report...", log)
			partitions2report_final = "all"

		# getting metadata report
		if args.metadata != "none":
			log.close()
			returned_value = run_metadata_report(args, partitions2report_final, "new")
			if str(returned_value) != "0":
				sys.exit("\n\nReporTree exited before expected while running metadata_report.py :-( Please check your input files and your command line. If you are in trouble and cannot figure out what is going on, contact us!!")
			cmds.append("metadata")
			log = open(log_name, "a+")
			
			# samples of interest
			s_of_interest = args.sample_of_interest
			if s_of_interest != "all":
				print_log("\tFiltering partitions_summary.tsv according to samples of interest...", log)
				samples_of_interest, singletons, do_not_exist = info_samples_interest(s_of_interest, args.output + "_partitions_summary.tsv", args.output + "_partitions.tsv", args.output)
		
		if "nomenclature" in cmds and "metadata" in cmds:
			nomenclature_change2summary(args.output)

	
	## others provided	--------------------
		
	elif args.allele_profile != "" or args.alignment != "" or args.vcf != "" or args.variants != "" or args.distance_matrix != "":
		distance_matrix_input = False
		if args.allele_profile != "": # ALLELE ---> DIRECT INPUT
			if args.alignment != "": 
				print_log("\nProfiles and sequence alignment files specified... I am confused :-(\n", log)
				sys.exit(1)
			
			elif args.vcf != "": 
				print_log("\nProfiles and VCF files specified... I am confused :-(\n", log)
				sys.exit(1)
			
			elif args.variants != "": 
				print_log("\nProfiles and variants list files specified... I am confused :-(\n", log)
				sys.exit(1)
			
			elif args.distance_matrix != "": 
				print_log("\nProfiles and distance matrix specified... I am confused :-(\n", log)
				sys.exit(1)
			
			else: # can continue using the allele profile
				analysis = args.analysis
				profile = args.allele_profile
				print_log("\nProfiles file provided -> will run partitioning_" + analysis + ".py:\n", log)
		
		elif args.alignment != "": # ALIGNMENT ---> PROCESS INPUT
			if args.vcf != "": 
				print_log("\nAlignment and VCF files specified... I am confused :-(\n", log)
				sys.exit(1)
			
			elif args.variants != "": 
				print_log("\nAlignment and variants list files specified... I am confused :-(\n", log)
				sys.exit(1)
			
			if args.distance_matrix != "": # needs to be elif
				print_log("\nAlignment and distance matrix specified... I am confused :-(\n", log)
				sys.exit(1)
			
			else: # can continue using the alignment
				analysis = args.analysis
				print_log("\nAlignment file provided -> will run alignment_processing.py and partitioning_" + analysis + ".py:\n", log)
				
				# processing alignment
				
				log.close()
				returned_value = run_alignment_processing(args)
				cmds.append("run_alignment_processing")
				if str(returned_value) != "0":
					sys.exit("\n\nReporTree exited before expected while running alignment_processing.py :-( Please check your input files and your command line. If you are in trouble and cannot figure out what is going on, contact us!!")
				log = open(log_name, "a+")
				if os.path.exists(args.output + "_align_profile.tsv"):
					profile = args.output + "_align_profile.tsv"
					input_align = True
				else:
					sys.exit(1)
				
		elif args.vcf != "": # VCF ---> PROCESS INPUT WITH VCF2MST
			if args.variants != "": 
				print_log("\nVCF and variants list files specified... I am confused :-(\n", log)
				sys.exit(1)
			elif args.distance_matrix != "":
				print_log("\nVCF and distance matrix specified... I am confused :-(\n", log)
				sys.exit(1)
			
			else: # can continue using the VCF
				analysis = args.analysis
				print_log("\nVCF file provided -> will run vcf2mst and partitioning_" + analysis + ".py:\n", log)
				print_log("\n-------------------- vcf2mst --------------------\n", log)
				print_log("perl " + reportree_path + "/scripts/vcf2mst/vcf2mst.pl " + args.vcf + " " + args.output + "_profile.tsv vcf -out profile", log)
				log.close()
				
				returned_value = run_vcf2mst(args, "vcf")
				cmds.append("vcf")
				if str(returned_value) != "0":
					sys.exit("\n\nReporTree exited before expected while running vcf2mst.pl :-( Please check your input files and your command line. If you are in trouble and cannot figure out what is going on, contact us!!")
				log = open(log_name, "a+")
				if os.path.exists(args.output + "_profile.tsv"):
					profile = args.output + "_profile.tsv"
				else:
					sys.exit(1)
					
		elif args.variants != "": # LIST ---> PROCESS INPUT WITH VCF2MST
			if args.distance_matrix != "":
				print_log("\nVariants list and distance matrix specified... I am confused :-(\n", log)
				sys.exit(1)
			else: # can continue using the variants list
				analysis = args.analysis
				print_log("\nVariants list provided -> will run vcf2mst and partitioning_" + analysis + ".py:\n", log)
				log.close()
				
				returned_value = run_vcf2mst(args, "var")
				cmds.append("var")
				if str(returned_value) != "0":
					sys.exit("\n\nReporTree exited before expected while running vcf2mst.pl :-( Please check your input files and your command line. If you are in trouble and cannot figure out what is going on, contact us!!")
				log = open(log_name, "a+")
				if os.path.exists(args.output + "_profile.tsv"):
					profile = args.output + "_profile.tsv"
				else:
					sys.exit(1)
		
		elif args.distance_matrix != "": # DISTANCE MATRIX ---> DIRECT INPUT
			analysis = "HC"
			distance_matrix_input = True
			profile = args.distance_matrix
			cmds.append("dist_mx")
			print_log("\nDistance matrix provided -> will run partitioning_" + analysis + ".py:\n", log)
			cmds.append(distance_matrix_input)
			
		
		# GRAPETREE	----------
		
		if analysis == "grapetree": # grapetree
			log.close()
			returned_value = run_partitioning_grapetree(args, input_align, profile)
			cmds.append("grapetree")
			if str(returned_value) != "0":
				sys.exit("\n\nReporTree exited before expected while running partitioning_grapetree.py :-( Please check your input files and your command line. If you are in trouble and cannot figure out what is going on, contact us!!")
			log = open(log_name, "a+")
		
		
		# HIERARCHICAL CLUSTERING	----------
		
		elif analysis == "HC": # hc	
			log.close()
			returned_value = run_partitioning_HC(args, input_align, distance_matrix_input, profile)
			cmds.append("HC")
			if str(returned_value) != "0":
				sys.exit("\n\nReporTree exited before expected while running partitioning_HC.py :-( Please check your input files and your command line. If you are in trouble and cannot figure out what is going on, contact us!!")
			log = open(log_name, "a+")
		
		else:
			print_log("\nAnalysis \"" + str(args.analysis) + "\" is not compatible with the provided input!!\n", log)
			sys.exit(1)
		

		# NOMENCLATURE	----------
		
		partition_status = "new"
		if args.nomenclature != "":
			print_log("\nNomenclature file provided. Will work on cluster names...", log)
			partitions = run_nomenclature(args.output + "_partitions.tsv", args.nomenclature, args.output, day, log)
			partition_status = "new"
			cmds.append("nomenclature")
		

		# running comparing partitions
		if "stability_regions" in args.partitions2report and "all" in args.partitions2report:
			print_log("\t'stability_regions' and 'all' options cannot be simultaneously specified in --partitions2report... I am confused :-(", log)
			sys.exit(1)
		elif "stability_regions" in args.partitions2report and "all" not in args.partitions2report:
			if args.threshold == "max" and args.dist == 1.0:
				print_log("\t'stability_regions' option specified. Will also run comparing_partitions_v2.py to determine the partitions that must be included in the report...", log)
				log.close()
				returned_value = run_stability(args, log_name, partition_status)
				cmds.append("stability")
				if str(returned_value) != "0":
					sys.exit("\n\nReporTree exited before expected while running comparing_partitions_v2.py :-( Please check your input files and your command line. If you are in trouble and cannot figure out what is going on, contact us!!")
			else:
				if args.threshold == "max" and args.dist != 1.0:
					print_log("\t'stability_regions' option specified but minimum distance was different from 1.", log)
					sys.exit(1)
				else:
					print_log("\t'stability_regions' option specified but partitions were not obtained for all possible thresholds.", log)
					sys.exit(1)
			log = open(log_name, "a+")
			
		partitions2report_final, partitions2report_lst = infer_partitions_from_list(analysis, args.output + "_partitions.tsv", metadata_col, args.partitions2report, args.output, args.dist, log)
		
		if len(partitions2report_final) == 0:
			print_log("\tThe analysis of the partitions to report returned an empty list. All partitions will be included in the report...", log)
			partitions2report_final = "all"
		
		# get nomenclature code
		if args.code_levels != "" and args.metadata != "":
			print_log("\nGenerating nomenclature code...", log)
			code_levels, code_levels_lst = infer_partitions_from_list(analysis, args.output + "_partitions.tsv", metadata_col, args.code_levels, args.output, args.dist, log)
			nomenclature_code(metadata_mx,metadata_col,args.output + "_partitions.tsv",args.output,code_levels_lst,day,log)
			args.metadata = args.output + "_metadata_w_partitions.tsv"

		# getting metadata report
		if args.metadata != "none":
			log.close()
			returned_value = run_metadata_report(args, partitions2report_final, "new")
			cmds.append("metadata_report")
			if str(returned_value) != "0":
				sys.exit("\n\nReporTree exited before expected while running metadata_report.py :-( Please check your input files and your command line. If you are in trouble and cannot figure out what is going on, contact us!!")
			log = open(log_name, "a+")

			# if allele matrix was filtered by loci called, add this info in metadata
			metadata = args.output + "_metadata_w_partitions.tsv"
			if os.path.exists(args.output + "_loci_report.tsv"):
				loci_called2metadata(metadata, args.output, args.loci_called, "loci")
			elif os.path.exists(args.output + "_pos_report.tsv"):
				loci_called2metadata(metadata, args.output, args.ATCG_content, "pos")

		# if nomenclature changed, add this info in partitions summary
		if "nomenclature" in cmds and "metadata_report" in cmds:
			nomenclature_change2summary(args.output)
			
		# samples of interest
		s_of_interest = args.sample_of_interest
		if s_of_interest != "all" and  args.metadata != "none":
			print_log("\n\n\n\n****************************** PROCESSING SAMPLES OF INTEREST ******************************\n\n", log)

			print_log("Filtering partitions_summary.tsv according to samples of interest...", log)
			samples_of_interest, singletons, do_not_exist = info_samples_interest(s_of_interest, args.output + "_partitions_summary.tsv", args.output + "_partitions.tsv", args.output)
			
			# UPDATE METADATA	----------
			metadata = interest2metadata(args.output, samples_of_interest)

			# ZOOM-IN	----------

			new_subset_filters = []

			if args.zoom != "no":
				print_log("Checking cluster zoom-in requests...", log)
				partitions4zoom, partitions4zoom_lst = infer_partitions_from_list(analysis, args.output + "_partitions.tsv", metadata_col, args.zoom, args.output, args.dist, log)
				clusters_of_interest, not_in_partitions = get_clusters_interest(samples_of_interest, args.output + "_partitions.tsv", args.output + "_metadata_w_partitions.tsv", day, partitions4zoom_lst)
				if len(partitions4zoom_lst) == 0:
					print_log("\tNone of the requested partitions for zoom-in is valid!... We are sorry but cluster zoom-in will be done! Please check that you have correctly indicated the list of partitions in the '--zoom-cluster-of-interest' argument.", log)
				else:
					for part in partitions4zoom_lst:
						if part in clusters_of_interest.keys():
							for cluster in clusters_of_interest[part]:
								tag_subset = str(part) + "_" + str(cluster)
								filter_subset = part + " == " + str(cluster)
								info = "zoom-in",tag_subset,filter_subset
								new_subset_filters.append(info)
						else:
							if part in not_in_partitions:
								print_log("\tClustering information was not obtained for the requested partition " + str(part) + ". So, we cannot do any zoom-in for this request. Please check that your '--zoom-cluster-of-interest' is coherent with your clustering requests.", log)
							else:
								print_log("\tAt the requested partition " + str(part) + ", we did not find any cluster with at least one sample of interest to zoom-in.", log)
			
			# TREE INTEREST	----------

			if args.subtree != "no" and args.metadata != "":
				print_log("Checking subtree requests...", log)
				n4subtree = args.subtree.split(",")
				hamming = args.output + "_dist_hamming.tsv"
				metadata = pandas.read_table(args.metadata, dtype = str)
				sample_col = metadata.columns[0]
				for n in n4subtree:
					for sample in samples_of_interest:
						closest_samples, code = get_closest_samples(sample, n, hamming)
						if code == "all_samples":
							print_log("\tSubtree of interest requested for n = " + str(n) + " will not be done because this includes all samples in the dataset!", log)
							break
						elif len(closest_samples) > 2:
							if code == "run":
								tag_subset = str(sample) + "_closest" + str(n)
								filter_subset = sample_col + " == " + ",".join(closest_samples)
								info = "subtree",tag_subset,filter_subset
								new_subset_filters.append(info)
							elif code == "not_in_hamming":
								print_log("\tSubtree of interest requested for n = " + str(n) + " for sample " + sample + " will not be done because this sample was not included in the initial analysis!", log)
			
			# SUBSET	----------

			zooms_file = open(args.output + "_zooms.txt", "w+")
			for type_analysis,tag_subset,filter_subset in new_subset_filters:
				out_folder = args.output
				if not os.path.exists(out_folder + "_" + tag_subset):
					os.system("mkdir " + out_folder + "_" + tag_subset)
				if "/" in out_folder:
					zooms_prefix = out_folder.split("/")[-1]
				else:
					zooms_prefix = out_folder
				print(zooms_prefix + "_" + tag_subset, file = zooms_file)
				input_align = False
				real_output =  args.output
				real_metadata = args.metadata
				real_filter = args.filter_column
				real_subset = args.subset
				args.metadata = metadata4subset(args.output + "_metadata_w_partitions.tsv", args.output + "_partitions.tsv", tag_subset, filter_subset)
				args.output = out_folder + "_" + tag_subset + "/" + tag_subset
				args.filter_column = filter_subset
				args.subset = True
				print_log("\n******************** WORKING ON: " + str(tag_subset) + " ********************", log)
				subset_log_name = args.output + ".log"
				subset_log = open(subset_log_name, "w+")
				print("\n******************** running reportree.py ********************\n", file = subset_log)
				print_log("version " + str(version) + " last updated on " + str(last_updated) + "\n", subset_log)
				print("\n-- THIS IS A SPECIAL RUN FOR SAMPLES OF INTEREST --\n", file = subset_log)
				print_log(" ".join(sys.argv), subset_log)
				subset_start = datetime.datetime.now()
				print_log("\nSubset start: " + str(start), subset_log)
				subset_log.close()
				cancel = False
				if "run_alignment_processing" in cmds:
					returned_value = run_alignment_processing(args)
					if str(returned_value) != "0":
						print_log(tag_subset + ": ReporTree cannot proceed after alignment_processing.py for this subset :-( Possibly you went out of samples or positions in the alignment.", log)
						cancel = True
					if os.path.exists(args.output + "_align_profile.tsv"):
						profile = args.output + "_align_profile.tsv"
						input_align = True
				if "grapetree" in cmds:
					if not cancel:
						args.threshold = "max"
						returned_value = run_partitioning_grapetree(args, input_align, profile)
						if str(returned_value) != "0":
							print_log(tag_subset + ": ReporTree cannot proceed after partitioning_grapetree.py :-( Possibly you went out of samples or positions/alleles in the matrix. Alternatively, this subset does not have diversity for partitions inference.", log)
							cancel = True
				elif "HC" in cmds:
					if not cancel:
						distance_matrix_input = False
						if "dist_mx" in cmds:
							distance_matrix_input = True
							profile = args.distance_matrix
						args.HCmethod_threshold = get_method_threshold(args.HCmethod_threshold)
						returned_value = run_partitioning_HC(args, input_align, distance_matrix_input, profile)
						if str(returned_value) != "0":
							print_log(tag_subset + ": ReporTree cannot proceed after partitioning_HC.py :-( Possibly you went out of samples or positions/alleles in the matrix. Alternatively, this subset does not have diversity for partitions inference.", log)
							cancel = True
				if "metadata_report" in cmds:
					if not cancel:
						rename_clusters_subsets(type_analysis,args.output + "_partitions.tsv",args.output)
						args.filter_column = ""
						args.subset = True
						returned_value = run_metadata_report(args, "all", "new")
						if str(returned_value) != "0":
							print_log(tag_subset + ": ReporTree cannot proceed after metadata_report.py :-( If you reached this step without previous errors, please contact us because this is a possible bug!!!", log)
				os.system("rm " + args.metadata)
				args.output = real_output
				args.metadata = real_metadata
				args.filter_column = real_filter
				args.subset = real_subset
				subset_log = open(subset_log_name, "a+")
				print_log("\n------------------------------------------------------------\n", subset_log)
				if not args.unzip:
					returned_value = os.system("zip -r " + out_folder + "_" + tag_subset + ".zip " + out_folder + "_" + tag_subset + "/")
					if str(returned_value) != "0":
						print_log("ReporTree failed to compress your cluster of interest or tree of interest directory: " + args.output + "_" + tag_subset + "/", subset_log)
					else:
						os.system("rm -rf " + out_folder + "_" + tag_subset + "/")
				subset_end = datetime.datetime.now()
				subset_elapsed = subset_end - subset_start
				print("ReporTree is done! If you found any issue please contact us!!\n", file = subset_log)
				print_log("\nSubset end: " + str(subset_end), subset_log)
				print_log("\nSubset time elapsed: " + str(subset_elapsed), subset_log)
				subset_log.close()
			zooms_file.close()

			
	## only metadata	--------------------
	
	elif args.metadata != "none": # only metadata provided
		print_log("\nOnly metadata file provided -> only metadata_report.py will be run:\n", log)
		log.close()
		partitions2report_final = args.partitions2report
		returned_value = run_metadata_report(args, partitions2report_final, "no")
		if str(returned_value) != "0":
			sys.exit("\n\nReporTree exited before expected :-( Please check your input files and your command line. If you are in trouble and cannot figure out what is going on, contact us!!")
		log = open(log_name, "a+")
	
	else:
		print("\nYOU NEED TO SPECIFY SOMETHING VALID... OTHERWISE, I DO NOT KNOW WHAT TO DO!!!!!!\n")
		sys.exit()
	
	
	# done	----------
	
	end = datetime.datetime.now()
	
	elapsed = end - start
	print_log("\n------------------------------------------------------------\n", log)
	print_log("ReporTree is done! If you found any issue please contact us!!\n", log)
	print_log("\nEnd: " + str(end), log)
	print_log("Time elapsed: " + str(elapsed), log)
	log.close()

if __name__ == "__main__":
    main()
