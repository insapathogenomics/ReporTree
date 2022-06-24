#!/usr/bin/env	python3

"""
Obtain genetic clusters at any partition level(s) of a distance matrixes using hierarchical clustering.

By Veronica Mixao
@INSA
"""

from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, maxdists, to_tree
from scipy.spatial.distance import squareform
import pandas
import argparse
import textwrap
import sys
import os


# functions	----------

def conv_nucl(alleles):
	"""convert nucl to integers"""
	
	mx = alleles

	alleles =  mx[mx.columns[1:]]
	alleles.insert(0, mx.columns[0], mx[mx.columns[0]])
	
	alleles.to_csv("temporary_profile.tsv", index = False, header=True, sep ="\t")
	
	
def filter_mx(matrix, mx, filters, matrix_type, log):
	""" filter the allele or pairwise distance matrix
	input: matrix
	output: filtered pandas dataframe
	"""
    
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
	
	if matrix_type == "dist":
		pairwise_dist_mx = pairwise_dist_mx.set_index(pairwise_dist_mx.columns[0], drop = True)
		
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
    max_dist = maxdists(clustering)[-1] + 1
    
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


def get_cluster_composition(partitions):
	""" get summary of cluster composition 
	input: clusters dataframe
	output: summary dataframe 
	"""
	
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
        
        
# running the pipeline	----------

if __name__ == "__main__":
    
	# argument options
    
	parser = argparse.ArgumentParser(prog="partitioning_HC.py", formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent("""\
									###############################################################################             
									#                                                                             #
									#                             partitioning_HC.py                              #
									#                                                                             #
									###############################################################################  
									                            
									partitioning_HC.py obtains genetic clusters at any distance threshold(s) of an
									allele/SNP profile or pairwise distance matrix using hierarchical clustering 
									methods.
									
									Note: if a profile is provided, pairwise hamming distances will be calculated.
									
									How to run partitioning_HC.py?
									
									A) Partitions at all thresholds using single linkage:
									partitioning_HC.py -d_mx DISTANCE_MATRIX -o OUTPUT_NAME --HC-threshold single 

									B) Partitions at single linkage threshold = 2 and average linkage
									in the range 10 to 30:
									partitioning_HC.py -d_mx DISTANCE_MATRIX -o OUTPUT_NAME --HC-threshold 
									single-2,average-10-30
									
									-----------------------------------------------------------------------------"""))
	
	group0 = parser.add_argument_group("Partitioning with Hierarchical Clustering", "Specifications to cut the tree with HC methods")
	group0.add_argument("-d_mx", "--distance_matrix", dest="distance_matrix", required=False, type=str, help="[OPTIONAL] Input pairwise distance matrix")
	group0.add_argument("-a", "--allele-profile", dest="allele_profile", required=False, type=str, help="[OPTIONAL] Input allele profile matrix (can either be an allele matrix or a SNP matrix)")
	group0.add_argument("-o", "--output", dest="out", required=True, type=str, help="[MANDATORY] Tag for output file name")
	group0.add_argument("--HC-threshold", dest="method_threshold", required=False, default="single", 
						help="List of HC methods and thresholds to include in the analysis (comma-separated). To get clustering at all possible thresholds for a given method, write \
						the method name (e.g. single). To get clustering at a specific threshold, indicate the threshold with a hyphen (e.g. single-10). To get clustering at a specific \
						range, indicate the range with a hyphen (e.g. single-2-10). Default: single (List of possible methods: single, complete, average, weighted, centroid, median, ward)")
	group0.add_argument("--loci-called", dest="loci_called", required=False, default = "", help="[OPTIONAL] Minimum percentage of loci called (e.g. '--loci-called 0.95' will only keep in the allele \
						matrix the samples with > 95%% of alleles called, i.e. <= 5%% missing data). Code for missing data: 0.")
	group0.add_argument("--site-inclusion", dest="samples_called", required=False, default = 0.0, help="[OPTIONAL: Useful to remove informative sites/loci with excess of missing data] Minimum \
						proportion of samples per site without missing data (e.g. '--site-inclusion 1.0' will only keep loci/positions without missing data, i.e. a core alignment; \
						'--site-inclusion 0.0' will keep all loci/positions) NOTE: This argument works on profile/alignment positions/loci (i.e. columns)! [default: 0.0]. Code for missing data: 0.")
	group0.add_argument("--wgMLST", dest="wgmlst", default=False, action="store_true", help="Set if your profile is based on wgMLST scheme (if set, '--loci-called' will be applied after \
						'--site-inclusion')")
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
	print(" ".join(sys.argv))
	print(" ".join(sys.argv), file = log)
	
	
	# processing allele profile ----------
	
	if args.allele_profile:
		print("Profile matrix provided... pairwise distance will be calculated!")
		print("Profile matrix provided... pairwise distance will be calculated!", file = log)
		
		allele_mx = pandas.read_table(args.allele_profile, dtype = str)
		allele_mx = allele_mx.replace({"N": "0", "a": "A", "c": "C", "t": "T", "g": "G"})
		
		
		# filtering allele matrix	----------
		
		if args.metadata != "" and args.filter_column != "":
			print("Filtering the distance matrix...")
			print("Filtering the distance matrix...", file = log)
			
			filters = args.filter_column
			mx = pandas.read_table(args.metadata, dtype = str)
			sample_column = mx.columns[0]
			initial_samples = len(allele_mx[allele_mx.columns[0]].values.tolist())
			allele_mx = filter_mx(allele_mx, mx, filters, "allele", log)
			final_samples = len(allele_mx[allele_mx.columns[0]].values.tolist())
			print("\tFrom the " + str(initial_samples) + " samples, " + str(final_samples) + " were kept in the matrix...")
			print("\tFrom the " + str(initial_samples) + " samples, " + str(final_samples) + " were kept in the matrix...", file = log)
			allele_mx.to_csv(args.out + "_subset_matrix.tsv", sep = "\t", index = None)
	
	
		# cleaning allele matrix (columns)	----------
		
		if float(args.samples_called) != 1.0:
			print("Keeping only sites/loci with information in >= " + str(float(args.samples_called) * 100) + "% of the samples...")
			print("Keeping only sites/loci with information in >= " + str(float(args.samples_called) * 100) + "% of the samples...", file = log)
			
			pos_t0 = len(allele_mx.columns[1:])
			for col in allele_mx.columns[1:]:
				values = allele_mx[col].values.tolist()
				if (len(values)-values.count("0"))/len(values) < float(args.samples_called):
					allele_mx = allele_mx.drop(columns=col)
			allele_mx.to_csv(args.out + "_flt_matrix.tsv", index = False, header=True, sep ="\t")
			pos_t1 = len(allele_mx.columns[1:])
			print("\tFrom the " + str(pos_t0) + " loci/positions, " + str(pos_t1) + " were kept in the matrix.")
			print("\tFrom the " + str(pos_t0) + " loci/positions, " + str(pos_t1) + " were kept in the matrix.", file = log)
		
		
		# cleaning allele matrix	----------

		if args.loci_called != "":
			print("Cleaning the profile matrix using a threshold of >" + str(args.loci_called) + " alleles/positions called per sample...")
			print("Cleaning the profile matrix using a threshold of >" + str(args.loci_called) + " alleles/positions called per sample...", file = log)
			
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
			allele_mx.to_csv(args.out + "_flt_matrix.tsv", index = False, header=True, sep ="\t")
			report_allele_df.to_csv(args.out + "_loci_report.tsv", index = False, header=True, sep ="\t")
			
			
		# getting distance matrix	----------
		
		print("Getting the pairwise distance matrix with cgmlst-dists...")
		print("Getting the pairwise distance matrix with cgmlst-dists...", file = log)
		
		
		# convert ATCG to integers
		conv_nucl(allele_mx)
		
		
		# run cgmlst-dists
		os.system("cgmlst-dists temporary_profile.tsv > " + args.out + "_dist.mx")
		os.system("rm temporary_profile.tsv")
		temp_df = pandas.read_table(args.out + "_dist.mx", dtype = str)
		temp_df.rename(columns = {"cgmlst-dists": "dists"}, inplace = True)
		temp_df.to_csv(args.out + "_dist.mx", sep = "\t", index = None)
		dist = pandas.read_table(args.out + "_dist.mx")
		
	
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
			sample_column = mx.columns[0]
			dist = filter_mx(dist, mx, filters, "dist", log)
			dist.to_csv(args.out + "_flt_dist.mx", sep = "\t", index = None)

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
	
	for combination in args.method_threshold.split(","):
		if "-" not in combination:
			method = combination
			threshold = "all"
		else:
			method = combination.split("-")[0]
			threshold = combination.split("-",1)[1]
		
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
		
		if threshold == "all":
			print("\tCalculating clustering in range",str(1),str(max_dist),"with a distance of",str(args.dist))
			print("\tCalculating clustering in range",str(1),str(max_dist),"with a distance of",str(args.dist), file = log)
			for thr in range(1,int(max_dist) + 1):
				partition = method + "-" + str(thr) + "x" + str(args.dist)
				clustering[partition] = list(fcluster(hc_matrix, t = int(thr) * args.dist, criterion = "distance"))
		else:
			if "-" in threshold:
				min_thr = int(threshold.split("-")[0])
				max_thr = int(threshold.split("-")[1]) + 1
				
				if max_thr > max_dist:
					max_thr = str(max_dist)
				
				print("\tCalculating clustering in range",str(min_thr),str(max_thr),"with a distance of",str(args.dist))
				print("\tCalculating clustering in range",str(min_thr),str(max_thr),"with a distance of",str(args.dist), file = log)
				for thr in range(min_thr,max_thr):
					partition = method + "-" + str(thr) + "x" + str(args.dist)
					clustering[partition] = list(fcluster(hc_matrix, t = thr * args.dist, criterion = "distance"))
			else:
				partition = method + "-" + str(threshold) + "x" + str(args.dist)
				print("\tCalculating clustering for threshold",str(threshold),"with a distance of",str(args.dist))
				print("\tCalculating clustering for threshold",str(threshold),"with a distance of",str(args.dist), file = log)
				clustering[partition] = list(fcluster(hc_matrix, t = int(threshold) * args.dist, criterion = "distance"))


	# output partitions
	
	print("Creating sample partitions file...")
	print("Creating sample partitions file...", file = log)
	df_clustering = pandas.DataFrame(data = clustering)
	df_clustering.to_csv(args.out + "_partitions.tsv", sep = "\t", index = None)
    
    
    # output cluster composition
	
	print("Creating cluster composition file...")
	print("Creating cluster composition file...", file = log)
	cluster_composition = get_cluster_composition(df_clustering)
	cluster_composition.to_csv(args.out + "_clusterComposition.tsv", index = False, header=True, sep ="\t")

print("\npartitioning_HC.py is done!")
print("\npartitioning_HC.py is done!", file = log)

log.close()
