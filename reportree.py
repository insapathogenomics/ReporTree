#!/usr/bin/env	python3

"""
Obtain partitions for minimum-spanning trees or rooted trees and a summary report of the clusters.

WARNING!! Currently, it is only able to deal with rooted trees or allele matrixes.

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
import glob


reportree_script = os.path.realpath(__file__)
reportree_path = reportree_script.rsplit("/", 1)[0]


# functions	----------

def get_partitions2report(analysis, partitions_field, output, distance):
	""" determines the columns of the partitions table
	to report. Depending on the type of analysis (treecluster,
	grapetree or other), this function interprets  the
	partitions2report argument.
	
	returns: 
	comma-separated list of columns for which a report
	is necessary 
	"""
	
	
	if partitions_field == "all":
		partitions2report = partitions_field
	else:
		partitions2report_lst = []
		if "," in partitions_field: # more than one method
			for part in partitions_field.split(","):
				if part == "stability_regions": # get stability regions
					for filename in glob.glob(output + "*stableRegions.tsv"):
						with open(filename, "r") as stable_blocks_open:
							stable_blocks = stable_blocks_open.readlines()
							for line in stable_blocks:
								if "#" not in line:
									l = line.split("\t")
									first_partition = l[1].split("->")[1]
									if first_partition not in partitions2report_lst:
										partitions2report_lst.append(first_partition)
						
				else:
					if analysis != "other": # run grapetree or treecluster
						if "-" in part: # threshold specified 
							if analysis == "grapetree":
								method_info = "MST"
								threshold_info = part
							else:
								method_info = part.split("-", 1)[0]
								threshold_info = part.split("-", 1)[1]
						
							if "-" in threshold_info: # range specified
								min_threshold = int(threshold_info.split("-")[0])
								max_threshold = int(threshold_info.split("-")[1]) + 1
						
								for threshold in range(min_threshold,max_threshold):
									info = method_info + "-" + str(threshold) + "x" + str(distance)
									partitions2report_lst.append(info)
							
							else: # specific number		
								info = method_info + "-" + str(threshold_info) + "x" + str(distance)
								partitions2report_lst.append(info)
						
						else: 
							if analysis == "grapetree": # specific number
								info = "MST-" + str(part) + "x" + str(distance)
								partitions2report_lst.append(info)
					else:
						partitions2report_lst.append(part)

		else: # single method
			if partitions_field == "stability_regions": # run grapetree or treecluster
				for filename in glob.glob(output + "*stableRegions.tsv"):
					with open(filename, "r") as stable_blocks_open:
						stable_blocks = stable_blocks_open.readlines()
						for line in stable_blocks:
							if "#" not in line:
								l = line.split("\t")
								first_partition = l[1].split("->")[1]
								if first_partition not in partitions2report_lst:
									partitions2report_lst.append(first_partition)
			else:
				if analysis != "other":
					if "-" in partitions_field: # threshold specified 
						if analysis == "grapetree":
							method_info = "MST"
							threshold_info = partitions_field
						else:
							method_info = partitions_field.split("-", 1)[0]
							threshold_info = partitions_field.split("-", 1)[1]
						
						if "-" in threshold_info: # range specified
							min_threshold = int(threshold_info.split("-")[0])
							max_threshold = int(threshold_info.split("-")[1]) + 1
						
							for threshold in range(min_threshold,max_threshold):
								info = method_info + "-" + str(threshold) + "x" + str(distance)
								partitions2report_lst.append(info)
							
						else: # specific number		
							info = method_info + "-" + str(threshold_info) + "x" + str(distance)
							partitions2report_lst.append(info)
						
					else: 
						if analysis == "grapetree": # specific number
							info = "MST-" + str(partitions_field) + "x" + str(distance)
							partitions2report_lst.append(info)
				else:
					partitions2report_lst.append(partitions_field)
	
		partitions2report = ",".join(partitions2report_lst)
	
	return partitions2report	


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
			

def filter_partitions_table(method, table):
	""" only keep the partitions table of the method in analysis """
	
	mx = pandas.read_table(table)
	
	suitable_columns = ["sequence"]
	for col in mx.columns.tolist():
		if col.split("-")[0] == method:
			suitable_columns.append(col)
	
	final_df = mx.filter(items=suitable_columns)
	
	with open("tmp.tsv", "w") as method_file:
		final_df.to_csv("tmp.tsv", index = False, header=True, sep ="\t")
		
		
def col_list(metadata, partitions, tree, allele_profile, method_threshold, root_dist_by_node, threshold, distance):
	""" output the list of columns 
	
	returns:
	list with the column names
	"""
	
	cols_output = []
	methods_output = []
	metadata_mx = pandas.read_table(metadata)
	
	for col in metadata_mx.columns:
		if col != "date":
			n_col = "n_" + col
			cols_output.append(col)
			cols_output.append(n_col)
		else:
			cols_output.append("first_seq_date")
			cols_output.append("last_seq_date")
			cols_output.append("timespan_days")
	
	if partitions != "":
		partitions_mx = pandas.read_table(partitions)
		for col in partitions_mx.columns:
			if col not in cols_output:
				n_col = "n_" + col
				cols_output.append(col)
				cols_output.append(n_col)
	
	elif tree != "":
		if "," in method_threshold:
			for mt in method_threshold.split(","):
				method = mt.split("-")[0]
				if method not in methods_output:
					methods_output.append(method)
		else:
			method = method_threshold.split("-")[0]
			methods_output.append(method)
			
	elif allele_profile != "":
		methods_output.append("MST")
	
	return cols_output,methods_output
	
	
# running the pipeline	----------

if __name__ == "__main__":
    
	# argument options
    
	parser = argparse.ArgumentParser(prog="reportree.py", formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent("""\
									#################################################################             
									#                                                               #
									#                         reportree.py                          #
									#                                                               #
									################################################################# 
									                            
									ReporTree aims to answer the need of:
									
									-obtaining genetic clusters at any partition level(s) of a rooted 
									tree (e.g. SNP-scaled tree) or of a minimum spanning tree derived 
									from cg/wgMLST allele matrix
									
									-obtaining summary reports with the statistics/trends (e.g. 
									timespan, location range, cluster/group size and composition, age 
									distribution etc.) for the derived genetic clusters or for any 
									other provided grouping variable (such as, clade, lineage, ST, 
									vaccination status, etc.)
									
									-identifying regions of cluster stability (i.e. cg/wgMLST 
									partition threshold ranges in which cluster composition is similar)
									
									
									Note: Column 'sequence' is mandatory in metadata and partition
									tables and white spaces should be avoided in column names!!
									
									Note 2: To take the most profit of this script we recommend that 
									you include the column 'date' in the metadata. This column must 
									follow the format YYYY-MM-DD. If you only provide YYYY, it will 
									assume YYYY-01-01!!
									
									Note 3: If a 'date' column is provided in the metadata, this 
									script will determine and provide in the new metadata table the 
									columns iso_year, iso_week_nr and iso_week for each sample  
									(e.g. iso_year = 2021, iso_week_nr = 52, iso_week = 2021W52)!!
									
									Note 4: While for nominal or categorical variables this script
									can provide in the summary report the number of observations or 
									the frequency of each observation, for the 'date' column this 
									script can provide:
										- first_seq_date
										- last_seq_date
										- timespan_days
									
									Note 5: Currently, for non-SNP-distance rooted trees, users 
									have to specify a minimum unit to cut the tree (currently, the
									default is 1, which is equivalent to 1 SNP in a SNP-scaled rooted
									tree). NEWS COMING SOON!!
									
									
									
									TIP!! If you do not know which columns you can indicate for the
									argument '--columns_summary_report', you can use the 
									'--list'!!
									
									
									NOTE: Check the github page for information about citation
									                  
									-----------------------------------------------------------------"""))
	
	
	## general parameters
	
	group0 = parser.add_argument_group("ReporTree", "ReporTree input/output file specifications")
	group0.add_argument("-m", "--metadata", dest="metadata", required=True, type=str, help="[MANDATORY] Metadata file in .tsv format (column 'sequence' is mandatory)")
	group0.add_argument("-t", "--tree", dest="tree", default="", required=False, type=str, help="[OPTIONAL] Input tree")
	group0.add_argument("-a", "--allele-profile", dest="allele_profile", default="", required=False, type=str, help="[OPTIONAL] Input allele profile matrix")
	group0.add_argument("-p", "--partitions", dest="partitions", required=False, default="", type=str, help="[OPTIONAL] Partitions file in .tsv format (column 'sequence' is mandatory) - \
						'partition' represents any variable that is not in the metadata")
	group0.add_argument("-out", "--output", dest="output", required=False, default="ReporTree", type=str, help="[OPTIONAL] Tag for output file name (default = ReporTree)")
	group0.add_argument("--list", dest="list_col_summary", required=False, action="store_true", help=" [OPTIONAL] If after your command line you specify this option, ReporTree will list all the \
						possible columns that you can use as input in '--columns_summary_report'. NOTE!! The objective of this argument is to help you with the input of '--columns_summary_report'. \
						So, it will not run reportree.py main functions!!")
	
	
	## partition threshold 
	
	group1 = parser.add_argument_group("Partition minimum unit", "Minimum unit/distance between partition thresholds")
	group1.add_argument("-d", "--dist", dest="dist", required=False, default=1.0, type=float, help="Distance unit by which partition thresholds will be multiplied (example: if -d 10 and \
						-thr 5,8,10-30, the minimum spanning tree will be cut at 50,80,100,110,120,...,300. If -d 10 and --method-threshold avg_clade-2, the avg_clade threshold will be set \
						at 20). This argument is particularly useful for non- SNP-scaled trees. Currently, the default is 1, which is equivalent to 1 allele distance or 1 SNP distance. [1.0]")
	
	## partitioning treecluster
	
	group1 = parser.add_argument_group("Partitioning with TreeCluster", "Specifications to cut the tree with TreeCluster [only if a tree file is provided]")
	group1.add_argument("--method-threshold", dest="method_threshold", required=False, default="root_dist,avg_clade-1", 
						help="List of TreeCluster methods and thresholds to include in the analysis (comma-separated). To get clustering at all possible thresholds for a given method, write \
						the method name (e.g. root_dist). To get clustering at a specific threshold, indicate the threshold with a hyphen (e.g. root_dist-10). To get clustering at a specific \
						range, indicate the range with a hyphen (e.g. root_dist-2-10). Default: root_dist,avg_clade-1 (List of possible methods: avg_clade, leaf_dist_max, leaf_dist_min, length, \
						length_clade, max, max_clade, root_dist, single_linkage, single_linkage_cut, single_linkage_union) Warning!! So far, ReporTree was only tested with avg_clade and \
						root_dist!")
	group1.add_argument("--support", dest="support", required=False, default="-inf", help="[OPTIONAL: see TreeCluster github for details] Branch support threshold") 
	group1.add_argument("--root-dist-by-node", dest="root_dist_by_node", required=False, action="store_false", help="[OPTIONAL] Set only if you WANT to cut the tree with root_dist method at each tree \
						node distance to the root (similar to root_dist at all levels but just for informative distances)")
	
	
	## partitioning grapetree
	
	group2 = parser.add_argument_group("Partitioning with GrapeTree", "Specifications to get and cut minimum spanning trees derived from cg/wgMLST allele data [only if an allele profile file is provided]")
	group2.add_argument("--method", dest="grapetree_method", default="MSTreeV2", help="\"MSTreeV2\" [DEFAULT]\n Alternative:\"MSTree\"\n")
	group2.add_argument("--missing", dest="handler", default=0, type=int, help="ONLY FOR MSTree. \n0: [DEFAULT] ignore missing data in pairwise comparison. \n1: remove column \
						with missing data. \n2: treat missing data as an allele. \n3: use absolute number of allelic differences.")
	group2.add_argument("--wgMLST", default=False, action="store_true", help="[EXPERIMENTAL: see GrapeTree github for details] a better support of wgMLST schemes")
	group2.add_argument("--n_proc",  dest="number_of_processes", type=int, default=5, help="Number of CPU processes in parallel use. [5]")
	group2.add_argument("-thr", "--threshold", dest="threshold", default = "max", help="Partition threshold for clustering definition. Different thresholds can be comma-separated (e.g. 5,8,16). \
						Ranges can be specified with a hyphen (e.g. 5,8,10-20). If this option is not set, the script will perform clustering for all the values in the range 1 to max")
	group2.add_argument("--subset", dest="subset", required=False, action="store_true", help="Reconstruct the minimum spanning tree using only the samples that correspond to the filters \
						specified at the '--filter' argument")	
	
	## reportree
	
	group3 = parser.add_argument_group("ReporTree metadata report", "Specific parameters to report clustering/grouping information associated to metadata")
	group3.add_argument("--columns_summary_report", dest="columns_summary_report", required=False, default="n_sequence,lineage,n_country,country,n_region,first_seq_date,last_seq_date,timespan_days", type=str, help="Columns \
						(i.e. variables of metadata) to get statistics for the derived genetic clusters or for other grouping variables defined in --metadata2report (comma-separated). If the \
						name of the column is provided, the different observations and the respective percentage are reported. If 'n_column' is specified, the number of the different \
						observations is reported. For example, if 'n_country' and 'country'  are specified, the summary will report the number of countries and their distribution (percentage) \
						per cluster/group. Exception: if a 'date' column is in the metadata, it can also report first_seq_date, last_seq_date, timespan_days. Default = \
						n_sequence,lineage,n_country,country,n_region,first_seq_date,last_seq_date,timespan_days [the order of the list will be the order of the columns in the report]")
	group3.add_argument("--partitions2report", dest="partitions2report", required=False, default="all", type=str, help="Columns of the partitions table to include in a joint report \
						(comma-separated). Other alternatives: 'all' == all partitions; 'stability_regions' == first partition of each stability region as determined by \
						comparing_partitions_v2.py. Warning!! 'stability_regions' can only be inferred when partitioning TreeCluster or GrapeTree is run for all possible thresholds or when a \
						similar partitions table is provided (i.e. sequential partitions obtained with the same clustering method) [all]")
	group3.add_argument("--metadata2report", dest="metadata2report", required=False, default="country", help="Columns of the metadata table for which a separated summary report must be \
						provided (comma-separated)")
	group3.add_argument("-f", "--filter", dest="filter_column", required=False, default="", help="[OPTIONAL] Filter for metadata columns to select the samples to analyze. This must be specified \
						within quotation marks in the following format 'column< >operation< >condition' (e.g. 'country == Portugal'). When more than one condition is specified for a given column, \
						they must be separated with commas (e.g 'country == Portugal,Spain,France'). When filters include more than one column, they must be separated with semicolon (e.g. \
						'country == Portugal,Spain,France;date > 2018-01-01;date < 2022-01-01'). White spaces are important in this argument, so, do not leave spaces before and after \
						commas/semicolons.")
	group3.add_argument("--frequency-matrix", dest="frequency_matrix", required=False, default="no", help="[OPTIONAL] Metadata column names for which a frequency matrix will be generated. This must \
						be specified within quotation marks in the following format 'variable1,variable2'. Variable1 is the variable for which frequencies will be calculated (e.g. for \
						'lineage,iso_week' the matrix reports the percentage of samples that correspond to each lineage per iso_week). If you want more than one matrix you can separate the \
						different requests with semicolon (e.g. 'lineage,iso_week;country,lineage'). If you want a higher detail in your variable2 and decompose it into two columns you use a colon \
						(e.g. lineage,country:iso_week will report the percentage of samples that correspond to each lineage per iso_week in each country)")
	group3.add_argument("--count-matrix", dest="count_matrix", required=False, default="no", help="[OPTIONAL] Same as '--frequency-matrix' but outputs counts and not frequencies")
	group3.add_argument("--mx-transpose", dest="mx_transpose", required=False, action="store_true", help="[OPTIONAL] Set ONLY if you want that the variable1 specified in '--frequency-matrix' \
						or in '--count-matrix' corresponds to the matrix first column.")
	group3.add_argument("--metadata-4-grapetree", dest="metadata4grapetree", required=False, action="store_true", help="Output an additional metadata file with the header ready for GrapeTree \
						visualization. Set only if you WANT the file")	
					
						
	## comparing partitions
	
	group4 = parser.add_argument_group("Stability regions", "Congruence analysis of cluster composition at all possible partitions to determine regions of cluster stability (automatically run if you set --partitions2report 'stability_regions'). WARNING! This option is planned to handle sequential partitions obtained with the same clustering method, such as a partitions table derived from cg/wgMLST data (from 1 to max allele threshold). Use it at your own risk, if you provide your own partitions table.")
	group4.add_argument("-AdjW", "--AdjustedWallace", dest="AdjustedWallace", action= "store", default=0.99, help="Threshold of Adjusted Wallace score to consider an observation for method stability \
						analysis [0.99]")
	group4.add_argument("-n", "--n_obs", dest="n_obs", action="store", default=5, help="Minimum number of sequencial observations that pass the Adjusted Wallace score to be considered a \
						'stability region' (i.e. a threshold range in which cluster composition is similar) [5]")
	group4.add_argument("-o", "--order", dest="order", action= "store", default=0, required=False, help="[Set only if you provide your own partitions table] Partitions order in the partitions \
						table (0: min -> max; 1: max -> min) [0]")
	
						 
	args = parser.parse_args()

	
	# check if the user wants the list of columns
	
	if args.list_col_summary:
		columns_metadata,columns_methods = col_list(args.metadata, args.partitions, args.tree, args.allele_profile, args.method_threshold, args.root_dist_by_node, args.threshold, args.dist)
		print("\nPossible column names to provide as input to '--columns_summary_report' argument:")
		print("\n".join(columns_metadata))
		
		if columns_methods != []:
			print("\nReporTree will also output some partition columns that can be included in your summary reports. For each method that you provided, the column name will have the following structure:")
			for method in columns_methods:
				print(method + "-<threshold>          ---------->          e.g. " + method + "-30")
				print(method + "-<range>          ---------->          e.g. " + method + "-30-40")
			print("\n")
		sys.exit()
		
		
	# starting logs	----------
    
	log_name = args.output + ".log"
	log = open(log_name, "w+")

	print("\n******************** running reportree.py ********************\n")
	print("\n******************** running reportree.py ********************\n", file = log)
	print(" ".join(sys.argv))
	print(" ".join(sys.argv), file = log)
	
	start = datetime.datetime.now()
	print("start:", start)
	print("start:", start, file = log)

    
    # reportree workflow	----------
    
    ## partitions table provided	--------------------
    
	if args.partitions != "": # partitions table provided
		if args.tree != "": # tree was provided -> need to get partitions again?
			print("\nTree and partitions files specified... I am confused :-(\n")
			print("\nTree and partitions files specified... I am confused :-(\n", file = log)
			sys.exit()
		elif args.allele_profile != "":
			print("\nAllele profiles and partitions files specified... I am confused :-(\n")
			print("\nAllele profiles and partitions files specified... I am confused :-(\n", file = log)
			sys.exit()
		else: # can continue using partitions
			print("\nPartitions file provided -> only metadata_report.py will be run:\n")
			print("\nPartitions file provided -> only metadata_report.py will be run:\n", file = log)
			
			# running comparing partitions
			if "stability_regions" in args.partitions2report and "all" in args.partitions2report:
				print("\t'stability_regions' and 'all' options cannot be simultaneously specified in --partitions2report... I am confused :-(")
				print("\t'stability_regions' and 'all' options cannot be simultaneously specified in --partitions2report... I am confused :-(", file = log)
				sys.exit()
			elif "stability_regions" in args.partitions2report and "all" not in args.partitions2report:
				print("\t'stability_regions' option specified. Will also run comparing_partitions_v2.py to determine the partitions that must be included in the report...")
				print("\t'stability_regions' option specified. Will also run comparing_partitions_v2.py to determine the partitions that must be included in the report...", file = log)
				log.close()
				os.system("python " + reportree_path + "/scripts/ComparingPartitions/comparing_partitions_v2.py -i1 " + args.partitions + " -t " + args.output + " -o1 " +  str(args.order) + " -a stability -n \
				" + str(args.n_obs) + " -thr " + str(args.AdjustedWallace) + " -log " + log_name)
				log = open(log_name, "a+")
				
			partitions2report_final = get_partitions2report("other", args.partitions2report, args.output, args.dist)
			
			if len(partitions2report_final) == 0:
				print("\tThe analysis of the partitions to report returned an empty list. All partitions will be included in the report...")
				print("\tThe analysis of the partitions to report returned an empty list. All partitions will be included in the report...", file = log)
				partitions2report_final = "all"
			log.close()
			
			# getting metadata report
			if args.metadata4grapetree:
				if args.mx_transpose:
					os.system("python " + reportree_path + "/scripts/metadata_report.py -m " + args.metadata + " -p " + args.partitions + " -o " + args.output + " --columns_summary_report \
					" + args.columns_summary_report + " --partitions2report " + partitions2report_final + " --metadata2report " + args.metadata2report + " -f \"" + args.filter_column + "\" \
					--metadata-4-grapetree --frequency-matrix \'" + args.frequency_matrix + "\' --count-matrix \'" + args.count_matrix + "\' --mx-transpose")
				else:
					os.system("python " + reportree_path + "/scripts/metadata_report.py -m " + args.metadata + " -p " + args.partitions + " -o " + args.output + " --columns_summary_report \
					" + args.columns_summary_report + " --partitions2report " + partitions2report_final + " --metadata2report " + args.metadata2report + " -f \"" + args.filter_column + "\" \
					--metadata-4-grapetree --frequency-matrix \'" + args.frequency_matrix + "\' --count-matrix \'" + args.count_matrix + "\'")
			else:
				if args.mx_transpose:
					os.system("python " + reportree_path + "/scripts/metadata_report.py -m " + args.metadata + " -p " + args.partitions + " -o " + args.output + " --columns_summary_report \
					" + args.columns_summary_report + " --partitions2report " + partitions2report_final + " --metadata2report " + args.metadata2report + " -f \"" + args.filter_column + "\" \
					--frequency-matrix \'" + args.frequency_matrix + "\' --count-matrix \'" + args.count_matrix + "\' --mx-transpose")
				else:
					os.system("python " + reportree_path + "/scripts/metadata_report.py -m " + args.metadata + " -p " + args.partitions + " -o " + args.output + " --columns_summary_report \
					" + args.columns_summary_report + " --partitions2report " + partitions2report_final + " --metadata2report " + args.metadata2report + " -f \"" + args.filter_column + "\" \
					--frequency-matrix \'" + args.frequency_matrix + "\' --count-matrix \'" + args.count_matrix + "\'")
			log = open(log_name, "a+")
	
	
	## tree provided	--------------------
	
	elif args.tree != "": # tree was provided
		if args.allele_profile != "": # allelic profiles provided -> grapetree or treecluster
			print("\nTree and allele profiles files specified... I am confused :-(\n")
			print("\nTree and allele profiles files specified... I am confused :-(\n", file = log)
			sys.exit()
		else: # can continue using tree and run treecluster and metadata report
			print("\nTree file provided -> will run partitioning_treecluster.py and metadata_report.py:\n")
			print("\nTree file provided -> will run partitioning_treecluster.py and metadata_report.py:\n", file = log)
			log.close()
			
			# running partitioning treecluster
			if args.root_dist_by_node == False:
				if args.support != "-inf":
					os.system("python " + reportree_path + "/scripts/partitioning_treecluster.py -t " + args.tree + " -o " + args.output + " --root-dist-by-node -d " + str(args.dist) + " \
					--method-threshold " + args.method_threshold + " --support " + str(args.support))
				else:
					os.system("python " + reportree_path + "/scripts/partitioning_treecluster.py -t " + args.tree + " -o " + args.output + " --root-dist-by-node -d " + str(args.dist) + " \
					--method-threshold " + args.method_threshold)
			else:
				if args.support != "-inf":
					os.system("python " + reportree_path + "/scripts/partitioning_treecluster.py -t " + args.tree + " -o " + args.output + " -d " + str(args.dist) + " --method-threshold \
					" + args.method_threshold + " --support " + str(args.support))
				else:
					os.system("python " + reportree_path + "/scripts/partitioning_treecluster.py -t " + args.tree + " -o " + args.output + " -d " + str(args.dist) + " --method-threshold \
					" + args.method_threshold)
			log = open(log_name, "a+")
		
		# running comparing partitions
		if "stability_regions" in args.partitions2report and "all" in args.partitions2report:
			print("\t'stability_regions' and 'all' options cannot be simultaneously specified in --partitions2report... I am confused :-(")
			print("\t'stability_regions' and 'all' options cannot be simultaneously specified in --partitions2report... I am confused :-(", file = log)
			sys.exit()
		elif "stability_regions" in args.partitions2report and "all" not in args.partitions2report:
			if len(all_partitions_available(args.method_threshold, args.root_dist_by_node)) >= 1 and args.dist == 1.0: # at least all partitions for 1 method
				print("\t'stability_regions' option specified. Will also run comparing_partitions_v2.py to determine the partitions that must be included in the report...")
				print("\t'stability_regions' option specified. Will also run comparing_partitions_v2.py to determine the partitions that must be included in the report...", file = log)
				for method in all_partitions_available(args.method_threshold, args.root_dist_by_node):
					filter_partitions_table(method, args.output + "_partitions.tsv")
					log.close()
					os.system("python " + reportree_path + "/scripts/ComparingPartitions/comparing_partitions_v2.py -i1 tmp.tsv -t " + args.output + "_" + method + " -o1 " + str(args.order) + " -a stability -n \
					" + str(args.n_obs) + " -thr " + str(args.AdjustedWallace) + " -log " + log_name)
					os.system("rm tmp.tsv")
			else:
				if len(all_partitions_available(args.method_threshold, args.root_dist_by_node)) >= 1 and args.dist != 1.0:
					print("\t'stability_regions' option specified but minimum distance was different from 1.")
					print("\t'stability_regions' option specified but minimum distance was different from 1.", file = log)
					sys.exit()
				else:
					print("\t'stability_regions' option specified but no method was run for all possible partitions.")
					print("\t'stability_regions' option specified but no method was run for all possible partitions.", file = log)
					sys.exit()
			
		partitions2report_final = get_partitions2report("treecluster", args.partitions2report, args.output, args.dist)
		
		if len(partitions2report_final) == 0:
			log = open(log_name, "a+")
			print("\tThe analysis of the partitions to report returned an empty list. All partitions will be included in the report...")
			print("\tThe analysis of the partitions to report returned an empty list. All partitions will be included in the report...", file = log)
			partitions2report_final = "all"
			log.close()

		# getting metadata report
		if args.mx_transpose:
			os.system("python " + reportree_path + "/scripts/metadata_report.py -m " + args.metadata + " -p " + args.output + "_partitions.tsv -o " + args.output + " --columns_summary_report \
			" + args.columns_summary_report + " --partitions2report " + partitions2report_final + " --metadata2report " + args.metadata2report + " -f \"" + args.filter_column + "\" \
			--frequency-matrix \'" + args.frequency_matrix + "\' --count-matrix \'" + args.count_matrix + "\' --mx-transpose")
		else:
			os.system("python " + reportree_path + "/scripts/metadata_report.py -m " + args.metadata + " -p " + args.output + "_partitions.tsv -o " + args.output + " --columns_summary_report \
			" + args.columns_summary_report + " --partitions2report " + partitions2report_final + " --metadata2report " + args.metadata2report + " -f \"" + args.filter_column + "\" \
			--frequency-matrix \'" + args.frequency_matrix + "\' --count-matrix \'" + args.count_matrix + "\'")
		log = open(log_name, "a+")


	## alleles provided	--------------------
		
	elif args.allele_profile != "": # allelic profiles provided
		print("\nAllelic profiles file provided -> will run partitioning_grapetree.py and metadata_report.py:\n")
		print("\nAllelic profiles file provided -> will run partitioning_grapetree.py and metadata_report.py:\n", file = log)
		log.close()
		
		# running partitioning grapetree
		if args.wgMLST == True:
			if args.subset == True:
				os.system("python " + reportree_path + "/scripts/partitioning_grapetree.py -a " + args.allele_profile + " -o " + args.output + " --method " + args.grapetree_method + " --missing \
				" + str(args.handler) + " --wgMLST --n_proc " + str(args.number_of_processes) + " -thr " + str(args.threshold) + " -d " + str(args.dist) + " -m " + args.metadata + " \
				-f \"" + args.filter_column + "\"")
			else:
				os.system("python " + reportree_path + "/scripts/partitioning_grapetree.py -a " + args.allele_profile + " -o " + args.output + " --method " + args.grapetree_method + " --missing \
				" + str(args.handler) + " --wgMLST --n_proc " + str(args.number_of_processes) + " -thr " + str(args.threshold) + " -d " + str(args.dist))
		else:
			if args.subset == True:
				os.system("python " + reportree_path + "/scripts/partitioning_grapetree.py -a " + args.allele_profile + " -o " + args.output + " --method " + args.grapetree_method + " --missing \
				" + str(args.handler) + " --n_proc " + str(args.number_of_processes) + " -thr " + str(args.threshold) + " -d " + str(args.dist) + " -m " + args.metadata + " \
				-f \"" + args.filter_column + "\"")
			else:
				os.system("python " + reportree_path + "/scripts/partitioning_grapetree.py -a " + args.allele_profile + " -o " + args.output + " --method " + args.grapetree_method + " --missing \
				" + str(args.handler) + " --n_proc " + str(args.number_of_processes) + " -thr " + str(args.threshold) + " -d " + str(args.dist))
		log = open(log_name, "a+")
		
		# running comparing partitions
		if "stability_regions" in args.partitions2report and "all" in args.partitions2report:
			print("\t'stability_regions' and 'all' options cannot be simultaneously specified in --partitions2report... I am confused :-(")
			print("\t'stability_regions' and 'all' options cannot be simultaneously specified in --partitions2report... I am confused :-(", file = log)
			sys.exit()
		elif "stability_regions" in args.partitions2report and "all" not in args.partitions2report:
			if args.threshold == "max" and args.dist == 1.0:
				print("\t'stability_regions' option specified. Will also run comparing_partitions_v2.py to determine the partitions that must be included in the report...")
				print("\t'stability_regions' option specified. Will also run comparing_partitions_v2.py to determine the partitions that must be included in the report...", file = log)
				log.close()
				os.system("python " + reportree_path + "/scripts/ComparingPartitions/comparing_partitions_v2.py -i1 " + args.output + "_partitions.tsv -t " + args.output + " -o1 " + str(args.order) + " -a \
				stability -n " + str(args.n_obs) + " -thr " + str(args.AdjustedWallace) + " -log " + log_name)
			else:
				if args.threshold == "max" and args.dist != 1.0:
					print("\t'stability_regions' option specified but minimum distance was different from 1.")
					print("\t'stability_regions' option specified but minimum distance was different from 1.", file = log)
					sys.exit()
				else:
					print("\t'stability_regions' option specified partitions were not obtained for all possible thresholds.")
					print("\t'stability_regions' option specified partitions were not obtained for all possible thresholds.", file = log)
					sys.exit()
			
		partitions2report_final = get_partitions2report("grapetree", args.partitions2report, args.output, args.dist)
				
		if len(partitions2report_final) == 0:
			log = open(log_name, "a+")
			print("\tThe analysis of the partitions to report returned an empty list. All partitions will be included in the report...")
			print("\tThe analysis of the partitions to report returned an empty list. All partitions will be included in the report...", file = log)
			partitions2report_final = "all"
			log.close()

		# getting metadata report
		if args.metadata4grapetree == True:
			if args.mx_transpose:
				os.system("python " + reportree_path + "/scripts/metadata_report.py -m " + args.metadata + " -p " + args.output + "_partitions.tsv -o " + args.output + " --columns_summary_report \
				" + args.columns_summary_report + " --partitions2report " + partitions2report_final + " --metadata2report " + args.metadata2report + " -f \"" + args.filter_column + "\" --metadata-4-grapetree\
				--frequency-matrix \'" + args.frequency_matrix + "\' --count-matrix \'" + args.count_matrix + "\' --mx-transpose")
			else:
				os.system("python " + reportree_path + "/scripts/metadata_report.py -m " + args.metadata + " -p " + args.output + "_partitions.tsv -o " + args.output + " --columns_summary_report \
				" + args.columns_summary_report + " --partitions2report " + partitions2report_final + " --metadata2report " + args.metadata2report + " -f \"" + args.filter_column + "\" --metadata-4-grapetree\
				--frequency-matrix \'" + args.frequency_matrix + "\' --count-matrix \'" + args.count_matrix + "\'")
		else:
			if args.mx_transpose:
				os.system("python " + reportree_path + "/scripts/metadata_report.py -m " + args.metadata + " -p " + args.output + "_partitions.tsv -o " + args.output + " --columns_summary_report \
				" + args.columns_summary_report + " --partitions2report " + partitions2report_final + " --metadata2report " + args.metadata2report + " -f \"" + args.filter_column + "\" \
				--frequency-matrix \'" + args.frequency_matrix + "\' --count-matrix \'" + args.count_matrix + "\' --mx-transpose")
			else:
				os.system("python " + reportree_path + "/scripts/metadata_report.py -m " + args.metadata + " -p " + args.output + "_partitions.tsv -o " + args.output + " --columns_summary_report \
				" + args.columns_summary_report + " --partitions2report " + partitions2report_final + " --metadata2report " + args.metadata2report + " -f \"" + args.filter_column + "\" \
				--frequency-matrix \'" + args.frequency_matrix + "\' --count-matrix \'" + args.count_matrix + "\'")
		log = open(log_name, "a+")

	
	## only metadata	--------------------
	
	else: # only metadata provided
		print("\nOnly metadata file provided -> only metadata_report.py will be run:\n")
		print("\nOnly metadata file provided -> only metadata_report.py will be run:\n", file = log)
		log.close()
		if args.metadata4grapetree == True:
			if args.mx_transpose:
				os.system("python " + reportree_path + "/scripts/metadata_report.py -m " + args.metadata + " -o " + args.output + " --columns_summary_report \
				" + args.columns_summary_report + " --partitions2report " + args.partitions2report + " --metadata2report " + args.metadata2report + " -f \"" + args.filter_column + "\" --metadata-4-grapetree \
				--frequency-matrix \'" + args.frequency_matrix + "\' --count-matrix \'" + args.count_matrix + "\' --mx-transpose")
			else:
				os.system("python " + reportree_path + "/scripts/metadata_report.py -m " + args.metadata + " -o " + args.output + " --columns_summary_report \
				" + args.columns_summary_report + " --partitions2report " + args.partitions2report + " --metadata2report " + args.metadata2report + " -f \"" + args.filter_column + "\" --metadata-4-grapetree \
				--frequency-matrix \'" + args.frequency_matrix + "\' --count-matrix \'" + args.count_matrix + "\'")
		else:
			if args.mx_transpose:
				os.system("python " + reportree_path + "/scripts/metadata_report.py -m " + args.metadata + " -o " + args.output + " --columns_summary_report \
				" + args.columns_summary_report + " --partitions2report " + args.partitions2report + " --metadata2report " + args.metadata2report + " -f \"" + args.filter_column + "\" \
				--frequency-matrix \'" + args.frequency_matrix + "\' --count-matrix \'" + args.count_matrix + "\' --mx-transpose")
			else:
				os.system("python " + reportree_path + "/scripts/metadata_report.py -m " + args.metadata + " -o " + args.output + " --columns_summary_report \
				" + args.columns_summary_report + " --partitions2report " + args.partitions2report + " --metadata2report " + args.metadata2report + " -f \"" + args.filter_column + "\" \
				--frequency-matrix \'" + args.frequency_matrix + "\' --count-matrix \'" + args.count_matrix + "\'")
		log = open(log_name, "a+")
	
	
	# done	----------
	
	end = datetime.datetime.now()
	
	elapsed = end - start
	print("\n------------------------------------------------------------\n")
	print("\n------------------------------------------------------------\n", file = log)
	print("ReporTree is done! If you found any issue please contact us!!\n")
	print("ReporTree is done! If you found any issue please contact us!!\n", file = log)
	print("\nEnd:", end)
	print("\nEnd:", end, file = log)
	print("Time elapsed:", elapsed)
	print("Time elapsed:", elapsed, file = log)
	log.close()
	
