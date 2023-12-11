#!/usr/bin/env	python3

"""
Obtain clustering information for multiple partitions of a newick tree using TreeCluster.py

Note: Currently, for non-SNP-distance trees, users have to provide a minimum distance!! NEWS COMING SOON!!
Note2: This script takes advantage of TreeCluster.py -> do not forget to cite its authors as well!!

By Veronica Mixao
@INSA
"""

import os
import sys
import argparse
import textwrap
import pandas
from datetime import date
import ete3 as ete

version = "1.3.0"
last_updated = "2023-12-11"

treecluster = "TreeCluster.py"


# functions	----------

def get_distances(nw, root):
	""" estimate distances """
	
	t = ete.Tree(nw)
	
	if root != "no":
		if root == "midpoint":
			root_point = t.get_midpoint_outgroup()
			t.set_outgroup(root_point)
		else:
			root_point = t.get_leaves_by_name(root)[0]
			t.set_outgroup(root_point)

	distances = []

	for node in t.traverse("levelorder"):
		dist = t.get_distance(node)
		if dist > 0:
			distances.append(float(dist))
	
	min_dist = min(distances)
	sort_dist = sorted(distances)
	farthest, max_dist = t.get_farthest_node()

	return float(min_dist), float(max_dist), set(sort_dist)
	

def parsing_method_threshold(include_node_partitions, method_threshold, dist, max_thr, distances, log):
	""" parsing method and threshold argument """
	
	list_runs = [] # method,threshold
	
	if include_node_partitions == True: # tree will be cut at each node
		distances = list(distances)
		for d in sorted(distances):
			info = "node",str(d)
			list_runs.append(info)
	
	if "," in method_threshold: # more than one method
		methods = method_threshold.split(",")
		for method in methods:
			if "-" in method: # threshold specified 
				method_info = method.split("-", 1)[0]
				threshold_info = method.split("-", 1)[1]
					
				if "-" in threshold_info: # range specified
					min_threshold = int(threshold_info.split("-")[0])
					max_threshold = int(threshold_info.split("-")[1]) + 1
					
					for threshold in range(min_threshold,max_threshold):
						if float(threshold) * float(dist) <= float(max_thr):
							info = method_info,str(threshold)
							list_runs.append(info)
						else:
							print("\t\t**Threshold " + str(threshold) + " for method " + str(method_info) + " will not be analyzed because it exceeds the maximum distance of the tree!**")
							print("\t\t**Threshold " + str(threshold) + " for method " + str(method_info) + " will not be analyzed because it exceeds the maximum distance of the tree!**", file = log)
				else: # specific number
					threshold = threshold_info			
					if float(threshold) * float(dist) <= float(max_thr):
						info = method_info,str(threshold_info)
						list_runs.append(info)
					else:
						print("\t\t**Threshold " + str(threshold) + " for method " + str(method_info) + " will not be analyzed because it exceeds the maximum distance of the tree!**")
						print("\t\t**Threshold " + str(threshold) + " for method " + str(method_info) + " will not be analyzed because it exceeds the maximum distance of the tree!**", file = log)
			else: # threshold not specified
				info = method,"all"
				list_runs.append(info)
	else: # single method
		if "-" in method_threshold: # threshold specified 
			method_info = method_threshold.split("-", 1)[0]
			threshold_info = method_threshold.split("-", 1)[1]
					
			if "-" in threshold_info: # range specified
				min_threshold = int(threshold_info.split("-")[0])
				max_threshold = int(threshold_info.split("-")[1]) + 1
					
				for threshold in range(min_threshold,max_threshold):
					if float(threshold) * float(dist)<= float(max_thr):
						info = method_info,str(threshold)
						list_runs.append(info)
					else:
						print("\t\t**Threshold " + str(threshold) + " for method " + str(method_info) + " will not be analyzed because it exceeds the maximum distance of the tree!**")
						print("\t\t**Threshold " + str(threshold) + " for method " + str(method_info) + " will not be analyzed because it exceeds the maximum distance of the tree!**", file = log)

			else: # range not specified
				threshold = threshold_info	
				if float(threshold) * float(dist) <= float(max_thr):
					info = method_info,str(threshold)
					list_runs.append(info)
				else:
					print("\t\t**Threshold " + str(threshold) + " for method " + str(method_info) + " will not be analyzed because it exceeds the maximum distance of the tree!**")
					print("\t\t**Threshold " + str(threshold) + " for method " + str(method_info) + " will not be analyzed because it exceeds the maximum distance of the tree!**", file = log)
		else: # threshold not specified
			info = method_threshold,"all"
			list_runs.append(info)
			
	return list_runs
	
	
def get_partitions(cluster_file, partitions, method, threshold, dist):
	""" create the partitions file and rename the -1 code by 'singleton' """
	
	mx = pandas.read_table(cluster_file, dtype = str)
	new_mx = mx.replace(str(-1), "singleton") # replacing singletons
	
	singleton_counter = 0
	for index, row in new_mx.iterrows():
		if row[new_mx.columns[1]] == "singleton":
			singleton_counter += 1
			row[new_mx.columns[1]] = "singleton_" + str(singleton_counter)
		else:
			row[new_mx.columns[1]] = "cluster_" + str(row[new_mx.columns[1]])
			
	if method == "node":
		column_name = method + "-" + str(threshold) + "x1.0" # defining column name
	else:
		column_name = method + "-" + str(threshold) + "x" + dist # defining column name
	new_mx_name = new_mx.rename(columns={"ClusterNumber": column_name})
	
	if partitions is not "":
		mx_partitions = pandas.DataFrame(data = partitions)
		a = mx_partitions.set_index("sequence", drop = False)
		b = new_mx_name.set_index("SequenceName")
		c = pandas.concat([a, b], axis=1)
		new_partitions = c.reset_index(drop = True)
	else:
		new_partitions = new_mx_name.rename(columns={"SequenceName": "sequence"})
	
	return new_partitions


def get_cluster_composition(partitions):
	""" get summary of cluster composition """
	
	summary = {"partition": [], "cluster": [], "cluster_length": [], "samples": []}
	order_columns = ["partition", "cluster", "cluster_length", "samples"]
	
	for col in partitions.columns:
		if col != "sequence":
			clusters = partitions[col].values.tolist()
			for cluster in set(clusters):
				flt_data = partitions[partitions[col] == cluster] # filter the dataframe
				summary["partition"].append(col)
				summary["cluster"].append(str(cluster))
				summary["cluster_length"].append(len(flt_data["sequence"].values.tolist()))
				summary["samples"].append(",".join(flt_data["sequence"].values.tolist()))
	
	summary_df = pandas.DataFrame(data = summary, columns = order_columns)
	
	return summary_df		

def filter_partitions_table(partitions, mx, filters, log):
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
	
	df = partitions[partitions[partitions.columns[0]].isin(samples)]
	
	return df	
	
# running the pipeline	----------

def main():
    
	# argument options
    
	parser = argparse.ArgumentParser(prog="partitioning_treecluster.py", formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent("""\
									###############################################################################             
									#                                                                             #
									#                         partitioning_treecluster.py                         #
									#                                                                             #
									############################################################################### 
									                            
									partitioning_treecluster.py obtains genetic clusters at any partition level(s) 
									of a newick tree (e.g. SNP-scaled tree)
									
									Note: Currently, for non-SNP-distance rooted trees, users have to specify a 
									minimum unit to cut the tree (the default is 1, which is equivalent to 1 SNP in 
									a SNP-scaled rooted tree). NEWS COMING SOON!!
									
									NOTE 2: Do not forget to also cite TreeCluster authors.
									
									-------------------------------------------------------------------------------"""))
	
	group0 = parser.add_argument_group("Partitioning with TreeCluster", "Specifications to cut the tree with TreeCluster")
	group0.add_argument("-t", "--tree", dest="tree", default="", required=True, type=str, help="[MANDATORY] Input newick tree")
	group0.add_argument("-o", "--output", dest="output", required=True, type=str, help="[MANDATORY] Tag for output file name")
	group0.add_argument("--method-threshold", dest="method_threshold", required=False, default="root_dist,avg_clade-1", 
						help="List of TreeCluster methods and thresholds to include in the analysis (comma-separated). To get clustering at all possible thresholds for a given method, write \
						the method name (e.g. root_dist). To get clustering at a specific threshold, indicate the threshold with a hyphen (e.g. root_dist-10). To get clustering at a specific \
						range, indicate the range with a hyphen (e.g. root_dist-2-10). Default: root_dist,avg_clade-1 (List of possible methods: avg_clade, leaf_dist_max, leaf_dist_min, length, \
						length_clade, max, max_clade, root_dist, single_linkage, single_linkage_cut, single_linkage_union) Warning!! So far, ReporTree was only tested with avg_clade and \
						root_dist!")
	group0.add_argument("--support", dest="support", required=False, default=float('-inf'), help="[OPTIONAL: see TreeCluster github for details] Branch support threshold") 
	group0.add_argument("--root-dist-by-node", dest="root_dist_by_node", required=False, action="store_true", help="[OPTIONAL] Set only if you WANT to cut the tree with root_dist method at \
						each tree node distance to the root (similar to root_dist at all levels but just for informative distances)")
	group0.add_argument("-d", "--dist", dest="dist", required=False, default=1.0, type=float, help="Distance unit by which partition thresholds will be multiplied (example: if -d 10 and \
						--method-threshold avg_clade-2, the avg_clade threshold will be set at 20). This argument is particularly useful for non- SNP-scaled trees to set the distance that may \
						correspond to 1 SNP difference. Currently, the default is 1, which is equivalent to 1 SNP distance. NEWS COMING SOON!![1.0]")
	group0.add_argument("-r", "--root", dest="root", required=False, default="no", help="Set root of the input tree. Specify the leaf name to use as output. Alternatively, write \
						'midpoint', if you want to apply midpoint rooting method.")
	group0.add_argument("-m", "--metadata", dest="metadata", required=False, default="", type=str, help="[OPTIONAL] Metadata file in .tsv format to generate the cluster composition file according \
		     			to the '--filter' argument")
	group0.add_argument("-f", "--filter", dest="filter_column", required=False, default="", help="[OPTIONAL] Filter for metadata columns to select the samples of the tree that must \
						be used to report cluster composition. This must be specified within quotation marks in the following format 'column< >operation< >condition' (e.g. 'country == Portugal'). When \
						more than one condition is specified for a given column, they must be separated with commas (e.g 'country == Portugal,Spain,France'). When filters include more than one \
						column, they must be separated with semicolon (e.g. 'country == Portugal,Spain,France;date > 2018-01-01;date < 2022-01-01'). White spaces are important in this argument, \
						so, do not leave spaces before and after commas/semicolons.")
	
	
	args = parser.parse_args()


	# starting logs	----------

	log_name = args.output + ".log"
	log = open(log_name, "a+")
	
	print("\n-------------------- partitioning_treecluster.py --------------------\n")
	print("\n-------------------- partitioning_treecluster.py --------------------\n", file = log)
	print("version", version, "last updated on", last_updated, "\n")
	print("version", version, "last updated on", last_updated, "\n", file = log)
	print(" ".join(sys.argv), "\n")
	print(" ".join(sys.argv), "\n", file = log)
	


	# preparing variables
	
	runs = []
	partitions = ""
	
	
	# assessing method_threshold info
	
	print("Assessing tree specificities...")
	print("Assessing tree specificities...", file = log)
	min_dist, max_dist, distances = get_distances(args.tree, args.root)
	min_dist = args.dist # default minimum distance is 1 because the script works better with snp trees
	
	print("\tThreshold magnitude is " + str(min_dist))
	print("\tThreshold magnitude is " + str(min_dist), file = log)
	print("\tMaximum possible threshold is " + str(max_dist))
	print("\tMaximum possible threshold is " + str(max_dist), file = log)
		
	list_runs = parsing_method_threshold(args.root_dist_by_node, args.method_threshold, min_dist, max_dist, distances, log) # parsing method_threshold input
	
	print("Getting partitions with the following command line(s):")
	print("Getting partitions with the following command line(s):", file = log)
	
	# TreeCluster loop
		
	for method,threshold in list_runs:
		if threshold != "all": # threshold is specific
			if method == "node":
				final_thr = threshold
				threshold = str(threshold)
				method_run = "root_dist"
				info_run = method + "-" + str(threshold) + "x1.0"
			else:
				final_thr = float(min_dist) * int(threshold)
				method_run = method
				info_run = method + "-" + str(threshold) + "x" + str(min_dist)
			cluster_file = args.output + ".tsv"
			print("\tTreeCluster.py -i " + args.tree + " -o " + cluster_file + " -t " + str(final_thr) + " -s " + str(args.support) + " -m " + method)
			print("\tTreeCluster.py -i " + args.tree + " -o " + cluster_file + " -t " + str(final_thr) + " -s " + str(args.support) + " -m " + method, file = log)
			
			if args.support != float('-inf'):
				returned_value = os.system("TreeCluster.py -i " + args.tree + " -o " + cluster_file + " -t " + str(final_thr) + " -s " + str(args.support) + " -m " + method_run)
				if str(returned_value) != "0":
					print("\nSomething went wrong while running TreeCluster :-( please double check your input files and ReporTree specifications!")
					print("\nSomething went wrong while running TreeCluster :-( please double check your input files and ReporTree specifications!", file = log)
					sys.exit(1)
			else:
				returned_value = os.system("TreeCluster.py -i " + args.tree + " -o " + cluster_file + " -t " + str(final_thr) + " -m " + method_run)
				if str(returned_value) != "0":
					print("\nSomething went wrong while running TreeCluster :-( please double check your input files and ReporTree specifications!")
					print("\nSomething went wrong while running TreeCluster :-( please double check your input files and ReporTree specifications!", file = log)
					sys.exit(1)
				
			partitions = get_partitions(cluster_file, partitions, method, threshold, str(min_dist))
			os.system("rm " + cluster_file)
			runs.append(info_run)
							
		else: # run for all possible thresholds
			i = 0
			while float(min_dist) * float(i) <= float(max_dist):
				final_thr = float(min_dist) * float(i)
				info_run = method + "-" + str(i) + "x" + str(min_dist)
				cluster_file = args.output + ".tsv"
				print("\tTreeCluster.py -i " + args.tree + " -o " + cluster_file + " -t " + str(final_thr) + " -s " + str(args.support) + " -m " + method)
				print("\tTreeCluster.py -i " + args.tree + " -o " + cluster_file + " -t " + str(final_thr) + " -s " + str(args.support) + " -m " + method, file = log)
				
				if args.support != float('-inf'):
					returned_value = os.system("TreeCluster.py -i " + args.tree + " -o " + cluster_file + " -t " + str(final_thr) + " -s " + str(args.support) + " -m " + method)
					if str(returned_value) != "0":
						print("\nSomething went wrong while running TreeCluster :-( please double check your input files and ReporTree specifications!")
						print("\nSomething went wrong while running TreeCluster :-( please double check your input files and ReporTree specifications!", file = log)
						sys.exit(1)
				else:
					returned_value = os.system("TreeCluster.py -i " + args.tree + " -o " + cluster_file + " -t " + str(final_thr) + " -m " + method)
					if str(returned_value) != "0":
						print("\nSomething went wrong while running TreeCluster :-( please double check your input files and ReporTree specifications!")
						print("\nSomething went wrong while running TreeCluster :-( please double check your input files and ReporTree specifications!", file = log)
						sys.exit(1)
				partitions = get_partitions(cluster_file, partitions, method, str(i), str(min_dist))
				os.system("rm " + cluster_file)
				runs.append(info_run)
							
				i += 1
							
				if float(min_dist) * float(i) > float(max_dist):
					print("\n\t***Stopping clustering for method " + method + " at threshold = " + str(final_thr) + " because the max threshold is " + str(max_dist) + "***\n")
					print("\n\t***Stopping clustering for method " + method + " at threshold = " + str(final_thr) + " because the max threshold is " + str(max_dist) + "***\n", file = log)


	# output partitions
	
	partitions = pandas.DataFrame(data = partitions)
	partitions.to_csv(args.output + "_partitions.tsv", index = False, header=True, sep ="\t")
	
	
	# output cluster composition
	if args.metadata != "" and args.filter_column != "":
		mx = pandas.read_table(args.metadata)
		partitions = filter_partitions_table(partitions, mx, args.filter_column, log)
	cluster_composition = get_cluster_composition(partitions)
	cluster_composition.to_csv(args.output + "_clusterComposition.tsv", index = False, header=True, sep ="\t")
	
	print("\npartitioning_treecluster.py is done!")
	print("\npartitioning_treecluster.py is done!", file = log)
	
	log.close()

if __name__ == "__main__":
    main()
