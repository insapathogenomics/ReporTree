from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, maxdists, to_tree

def dist_mx(dist, logger):
    """ obtain a condensed distance matrix
    input: squared pairwise distance matrix
    output: consensed distance matrix
    """
    
    logger.info("Getting condensed distance matrix...")

    dist = dist.set_index(dist.columns[0],drop=True)
    
    samples = dist.columns
    
    condensed_dist = squareform(dist)
    
    return dist, condensed_dist, samples

    			
def hcluster(dist_mx, method_choice, logger):
    """ obtain linkage array 
    input: distance matrix
    output: linkage array, list samples names, max dist
    """
    
    clustering = linkage(dist_mx, method = method_choice)
    max_dist = maxdists(clustering)[-1]
    
    logger.info("Maximum distance " + str(max_dist) + "...")
    
    return clustering, max_dist
	

def get_partitions(clustering, threshold, logger):    
    """ obtain clustering information for a given threshold
    input: linkage array
    output: list with cluster names
    """
    
    clusters = list(fcluster(clustering, t = threshold, criterion = "distance"))
    
    return clusters


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


def hierarchical_clustering(df_dist, logger, args):

	distance_mx, condensed_dist_mx, samples = dist_mx(df_dist, logger)
	
	clustering = {"sequence": distance_mx.columns.tolist()}
	
	combinations2run = {}
	pct_correspondence = {}
	
	if args.pct_HCmethod_threshold != "none":
		logger.info("Correspondence between percentage and number of differences:")
		logger.info("#METHODPERCENTAGEDIFFERENCES")
		for combination_pct in args.pct_HCmethod_threshold.split(","):
			method = combination_pct.split("-")[0]
			threshold_pct = combination_pct.split("-",1)[1]
			threshold = str(int(int(total_size) * float(threshold_pct)))
			info_run = threshold,"pct"
			
			if method not in combinations2run.keys():
				combinations2run[method] = []
			combinations2run[method].append(info_run)

			if threshold not in pct_correspondence.keys():
				pct_correspondence[threshold] = []
			if str(threshold_pct) not in pct_correspondence[threshold]:
				pct_correspondence[threshold].append(str(threshold_pct))
				logger.info("" + str(float(threshold_pct)*100) + "" + str(threshold))
		
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
		logger.info("Hierarchical clustering with method: " + method + "...")
		hc_matrix, max_dist = hcluster(condensed_dist_mx, method, logger)
		
		# get newick
		
		logger.info("Generating newick file...")
	
		tree = to_tree(hc_matrix, False)
		nw = get_newick(tree, tree.dist, samples)
		
		with open(args.out + "_" + method + "_HC.nwk", "w+") as newick_out:
			print(nw, file = newick_out)
		
		# partitioning
		
		logger.info("Defining clusters...")
		
		
		for threshold,request in combinations2run[method]:
			if threshold == "all":
				logger.info(f"Calculating clustering in range 0.0 - {str(max_dist)} with a distance of {str(args.dist)}")
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
					
					if max_thr > max_dist:
						max_thr = str(max_dist)
					
					logger.info("Calculating clustering in range",str(min_thr),str(max_thr),"with a distance of",str(args.dist))
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
						partition = method + "-" + str(threshold) + "x" + str(args.dist)
						if partition not in cluster_details.keys():
							cluster_details[partition] = {}
						logger.info("Calculating clustering for threshold",str(threshold),"with a distance of",str(args.dist))
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
						partition = method + "-" + str(threshold) + "_(" + "_".join(pct_correspondence[threshold]) + ")"
						if partition not in cluster_details.keys():
							cluster_details[partition] = {}
						logger.info("Calculating clustering for threshold " + method + "-" + str(threshold) + ", which corresponds to the pct threshold of: " + ", ".join(pct_correspondence[threshold]))
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
	return clustering, cluster_details
