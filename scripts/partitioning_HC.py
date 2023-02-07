#!/usr/bin/env	python3

"""
Obtain genetic clusters at any partition level(s) of a distance matrixes using hierarchical clustering
By Veronica Mixao
@INSA

Modified by Finn Gruwier Larsen <figl@ssi.dk>
"""

import sys
import os
import argparse
import textwrap
import pandas
from datetime import date
import logging
import subprocess
from io import StringIO
from pathlib import Path
import uuid

from hierarchical_clustering import hierarchical_clustering

version = "1.1.2_ssi"
last_updated = "2023-02"

def create_logger(out: str):
	fh = logging.FileHandler(out + '.log', mode='a')
	fh.setLevel(logging.DEBUG)
	ch = logging.StreamHandler(sys.stdout)
	ch.setLevel(logging.INFO)
	logger = logging.getLogger()
	logger.setLevel(logging.INFO)
	logger.addHandler(fh)
	logger.addHandler(ch)
	return logger

def generate_random_filename():
    return uuid.uuid4().hex[:8].upper()

class HC:
	"""
	Instantiate this class when importing partitioning_HC into another code module;
	then use the 'run' method to run the calculation.
	"""
	out:str
	distance_matrix:str = ''
	allele_profile:str = ''
	method_threshold:str = 'single'
	pct_HCmethod_threshold: str = 'none'
	samples_called:float = 0.0
	loci_called:float = ''
	metadata:str = ''
	filter_column:str = ''
	dist:float = 1.0
	df_dist:pandas.DataFrame = None

	def __init__(self, out, **kwargs):
		self.out = out
		self.__dict__.update(kwargs)
	
	def __str__(self):
		return f"Output file prefix: {self.out}\n" + \
			f"Distance matrix filename: {self.distance_matrix}\n" + \
			f"Allele profile filename: {self.allele_profile}\n" + \
			f"Metadata filename: {self.metadata}\n" + \
			f"Method threshold: {self.method_threshold}\n" + \
			f"Percentage method threshold: {self.pct_HCmethod_threshold}\n" + \
			f"Samples called: {self.samples_called}\n" + \
			f"Loci called: {self.loci_called}\n" + \
			f"Filter column: {self.filter_column}\n" + \
			f"Distances: {self.dist}"
	
	def run(self):
		self.logger = create_logger(self.out)
		self.logger.info("Running hierarchical clustering with these parameters:")
		self.logger.info(self.__str__())
		self.logger.info("")

		if self.allele_profile:
			self.logger.info("Profile matrix provided; pairwise distance will be calculated.")
			self.df_dist = from_allele_profile(self, self.logger)
		elif self.distance_matrix:
			self.logger.info("Distance matrix provided; pairwise distance will NOT be calculated.")
			self.df_dist = from_distance_matrix(self, self.logger)
		else:
			msg = "Either distance matrix or allele profile must be specified!"
			self.logger.error(msg)
			raise ValueError(msg)
		
		# Copy-pasted from main part
		clustering, cluster_details = hierarchical_clustering(self.df_dist, self.logger, self)

		# output partitions
	
		self.logger.info("Creating sample partitions file...")
		df_clustering = pandas.DataFrame(data = clustering)
		df_clustering.to_csv(args.out + "_partitions.tsv", sep = "\t", index = None)
		
		
		# output cluster composition
		
		self.logger.info("Creating cluster composition file...")
		cluster_composition = get_cluster_composition(args.out + "_clusterComposition.tsv", cluster_details)
		#cluster_composition.to_csv(args.out + "_clusterComposition.tsv", index = False, header=True, sep ="\t")

		self.logger.info("partitioning_HC.py is done!")


# functions	----------

def conv_nucl(alleles):
	"""convert nucl to integers"""
	
	alleles = alleles.replace({"N": "0", "A": "1", "C": "2", "T": "3", "G": "4"})
	return alleles
	
	
def filter_mx(matrix, mx, filters, matrix_type, logger):
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
		mx["iso_week_nr"] = week.astype(str)
		mx["iso_week"] = year.astype(str).replace("<NA>", "-") + "W" + week.astype(str).replace("<NA>", "--")
		isoyear = mx.pop("iso_year")
		isoweek = mx.pop("iso_week_nr")
		isodate = mx.pop("iso_week")
		mx.insert(index_no + 1, "iso_year", isoyear)
		mx.insert(index_no + 2, "iso_week_nr", isoweek)
		mx.insert(index_no + 3, "iso_week", isodate)
				
	logger.info("\tFiltering metadata for the following parameters: " + " & ".join(filters.split(";")))
			
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

def from_allele_profile(hc=None, logger=None):
		global args
		if hc:
			args = hc
		allele_mx = pandas.read_table(args.allele_profile, dtype = str)
		allele_mx = allele_mx.replace("INF-","", regex=True) #implemented because of chewie-ns profiles
		allele_mx = allele_mx.replace("\*","", regex=True) #implemented because of chewie-ns profiles
		allele_mx = allele_mx.replace({"N": "0", "a": "A", "c": "C", "t": "T", "g": "G"})
		
		
		# filtering allele matrix	----------
		
		if args.metadata and args.filter_column:
			logger.info("Filtering the distance matrix...")
			
			filters = args.filter_column
			mx = pandas.read_table(args.metadata, dtype = str)
			sample_column = mx.columns[0]
			initial_samples = len(allele_mx[allele_mx.columns[0]].values.tolist())
			allele_mx = filter_mx(allele_mx, mx, filters, "allele", logger)
			final_samples = len(allele_mx[allele_mx.columns[0]].values.tolist())
			logger.info("\tFrom the " + str(initial_samples) + " samples, " + str(final_samples) + " were kept in the matrix...")
			allele_mx.to_csv(args.out + "_subset_matrix.tsv", sep = "\t", index = None)
	
	
		# cleaning allele matrix (columns)	----------
		
		if args.samples_called and float(args.samples_called) != 1.0 and float(args.samples_called) != 0.0:
			logger.info("Keeping only sites/loci with information in >= " + str(float(args.samples_called) * 100) + "% of the samples...")
			
			pos_t0 = len(allele_mx.columns[1:])
			for col in allele_mx.columns[1:]:
				values = allele_mx[col].values.tolist()
				if (len(values)-values.count("0"))/len(values) < float(args.samples_called):
					allele_mx = allele_mx.drop(columns=col)
			allele_mx.to_csv(args.out + "_flt_matrix.tsv", index = False, header=True, sep ="\t")
			pos_t1 = len(allele_mx.columns[1:])
			logger.info("\tFrom the " + str(pos_t0) + " loci/positions, " + str(pos_t1) + " were kept in the matrix.")
		
		
		# cleaning allele matrix (rows)	----------

		if args.loci_called:
			logger.info("Cleaning the profile matrix using a threshold of >" + str(args.loci_called) + " alleles/positions called per sample...")
			
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
			
			if len(pass_samples) == 0:
				logger.info("Cannot proceed because " + str(len(pass_samples)) + " samples were kept in the matrix!")
				sys.exit()
		
			logger.info("\tFrom the " + str(len(allele_mx[allele_mx.columns[0]].values.tolist())) + " samples, " + str(len(pass_samples)) + " were kept in the profile matrix.")
			
			allele_mx = allele_mx[allele_mx[allele_mx.columns[0]].isin(pass_samples)]
			allele_mx.to_csv(args.out + "_flt_matrix.tsv", index = False, header=True, sep ="\t")
			report_allele_df.to_csv(args.out + "_loci_report.tsv", index = False, header=True, sep ="\t")
			
			
		# getting distance matrix	----------
		
		logger.info("Getting the pairwise distance matrix with cgmlst-dists...")
		
		
		# convert ATCG to integers
		allele_mx = conv_nucl(allele_mx)


		# save allele matrix to a file that cgmlst-dists can use for input
		tmp_dirname = os.getenv('TMPDIR', '/tmp')
		tmp_path = Path(tmp_dirname, generate_random_filename())
		allele_mx.to_csv(tmp_path, index = False, header=True, sep ="\t")
		total_size = len(allele_mx.columns) - 1
		
		
		# run cgmlst-dists
		cp1:subprocess.CompletedProcess = subprocess.run(
			["cgmlst-dists", str(tmp_path)], capture_output=True, text=True)
		if cp1.returncode != 0:
			logger.error(f"Could not run cgmlst-dists on {str(tmp_path)}!")
			logger.error(cp1.stderr)
			raise IOError
		
		# remove temporary file
		cp2:subprocess.CompletedProcess = subprocess.run(
			["rm", str(tmp_path)], capture_output=True, text=True)
		logger.info(cp2.stdout)
		if cp2.returncode != 0:
			logger.warning(f"Could not remove temporary file {str(tmp_path)}!")
			logger.warning(cp2.stderr)
		
		temp_df = pandas.read_table(StringIO(cp1.stdout), dtype=str)
		temp_df.rename(columns = {"cgmlst-dists": "dists"}, inplace = True)
		temp_df.to_csv(args.out + "_dist.tsv", sep = "\t", index = None)
		dist = pandas.read_table(args.out + "_dist.tsv")
		return dist

def from_distance_matrix(hc=None, logger=None):
	global args
	if hc:
		args = hc
	dist = pandas.read_table(args.distance_matrix)

	# filtering the pairwise distance matrix	----------
	
	if args.metadata and args.filter_column:
		logger.info("Filtering the distance matrix...")
		
		filters = args.filter_column
		mx = pandas.read_table(args.metadata, dtype = str)
		sample_column = mx.columns[0]
		dist = filter_mx(dist, mx, filters, "dist", logger)
		dist.to_csv(args.out + "_flt_dist.tsv", sep = "\t", index = None)

	elif args.metadata and args.filter_column == "":
		logger.info("Metadata file was provided but no filter was found... I am confused :-(")
		sys.exit()

	elif (not args.metadata) and args.filter_column:
		logger.info("Metadata file was not provided but a filter was found... I am confused :-(")
		sys.exit()

	else:
		sample_column = "sequence"

	return dist


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
									methods. If a profile is provided, pairwise hamming distances will be 
									calculated.
									
									Note: the profile matrix provided can only contain integers or ATCG values.
									Alleles starting by "*" or "INF-" will be considered as a new. All the other 
									values will be considered as missing data.
									
									-----------------------------------------------------------------------------"""))
	
	group0 = parser.add_argument_group("Partitioning with Hierarchical Clustering", "Specifications to get partitions with HC methods")
	group0.add_argument("-d_mx", "--distance_matrix", dest="distance_matrix", required=False, type=str, help="[OPTIONAL] Input pairwise distance matrix")
	group0.add_argument("-a", "--allele-profile", dest="allele_profile", required=False, type=str, help="[OPTIONAL] Input allele profile matrix (can either be an allele matrix or a SNP matrix)")
	group0.add_argument("-o", "--output", dest="out", required=True, type=str, help="[MANDATORY] Tag for output file name")
	group0.add_argument("--HC-threshold", dest="method_threshold", required=False, default="single", 
						help="[OPTIONAL] List of HC methods and thresholds to include in the analysis (comma-separated). To get clustering at all possible thresholds for a given method, just write \
						the method name (e.g. single). To get clustering at a specific threshold, indicate the threshold with a hyphen (e.g. single-10). To get clustering at a specific \
						range, indicate the range with a hyphen (e.g. single-2-10). Note: Threshold values are inclusive, i.e. '--HC-threshold single-7' will consider samples with <= 7 \
						differences as belonging to the same cluster! Default: single (List of possible methods: single, complete, average, weighted, centroid, median, ward)")
	group0.add_argument("--pct-HC-threshold", dest="pct_HCmethod_threshold", required=False, default="none", help="[OPTIONAL] Similar to '--HC-threshold' but the partition threshold for cluster definition \
						is set as the proportion of differences to the final allelic schema size or number of informative positions, e.g. '--pct-HC-threshold single-0.005' corresponds to a \
						threshold of 5 allelic/SNP differences in a matrix with 1000 loci/sites under analysis. Ranges CANNOT be specified.")
	group0.add_argument("--site-inclusion", dest="samples_called", required=False, default = 0.0, help="[OPTIONAL: Useful to remove informative sites/loci with excess of missing data] Minimum \
						proportion of samples per site/loci without missing data (e.g. '--site-inclusion 1.0' will only keep loci/positions without missing data, i.e. a core alignment; \
						'--site-inclusion 0.0' will keep all loci/positions) NOTE: This argument works on profile/alignment loci/positions (i.e. columns)! [default: 1.0]. Code for missing data: 0.")
	group0.add_argument("--loci-called", dest="loci_called", required=False, default = "", help="[OPTIONAL] Minimum percentage of loci/positions called for allele/SNP matrices (e.g. \
						'--loci-called 0.95' will only keep in the profile matrix samples with > 95%% of alleles/positions, i.e. <= 5%% missing data). Applied after '--site-inclusion' argument! \
						Code for missing data: 0.")
	group0.add_argument("-m", "--metadata", dest="metadata", required=False, default="", type=str, help="[OPTIONAL] Metadata file in .tsv format to select the samples to use for clustering \
						according to the '--filter' argument")
	group0.add_argument("-f", "--filter", dest="filter_column", required=False, default="", help="[OPTIONAL] Filter for metadata columns to select the samples that must be used for HC \
						clustering. This must be specified within quotation marks in the following format 'column< >operation< >condition' (e.g. 'country == Portugal'). When \
						more than one condition is specified for a given column, they must be separated with commas (e.g 'country == Portugal,Spain,France'). When filters include more than one \
						column, they must be separated with semicolon (e.g. 'country == Portugal,Spain,France;date > 2018-01-01;date < 2022-01-01'). White spaces are important in this argument, \
						so, do not leave spaces before and after commas/semicolons.")
	group0.add_argument("-d", "--dist", dest="dist", required=False, default=1.0, type=float, help="Distance unit by which partition thresholds will be multiplied (example: if -d 10 and \
						--method-threshold single-2, the single linkage threshold will be set at 20).")
	
	
	hc = HC(**vars(parser.parse_args()))
	hc.run()