"""
Runs Several Tests for ReporTree functionality
author daniel.sobral@insa.min-saude.pt
"""
import os
import shutil
import subprocess
import pandas as pd
import datetime

def test_alignment():
	
	if not os.path.exists(os.path.join("tests","TEST1")):
		os.mkdir(os.path.join("tests","TEST1"))

	result = subprocess.check_output("python reportree.py " + \
		"-align tests/alignment.fas -out tests/TEST1/TEST1 --analysis grapetree -m tests/metadata.tsv " + \
		"-f \'country_origin != B;date >= 2023-02-02\' --subset --missing-code N --site-inclusion 0.9 --sample-ATCG-content  0.9 " + \
		"--remove-reference --use-reference-coords -r sample_reference --get-position-correspondence all " + \
		"--nomenclature-file tests/nomenclature.tsv --nomenclature-code-levels 15,7,4 --columns_summary_report " + \
		"first_seq_date,last_seq_date,timespan_days,country_origin,n_country_origin,source_isolation,n_source_isolation " + \
		"--partitions2report stability_regions --sample_of_interest sample9 --zoom-cluster-of-interest 15 --subtree-of-interest 5 " + \
		"--frequency-matrix country_origin,source_isolation --count-matrix nomenclature_code,country_origin --pivot", shell=True)    

    

    ### Check if files exist:
	assert(os.path.exists(os.path.join("tests","TEST1","TEST1.log")))
	assert(os.path.exists(os.path.join("tests","TEST1","TEST1.nwk")))
	assert(os.path.exists(os.path.join("tests","TEST1","TEST1_SAMPLES_OF_INTEREST_partitions_summary.tsv")))
	assert(os.path.exists(os.path.join("tests","TEST1","TEST1_MST-15x1.0_cluster_2.zip")))
	assert(os.path.exists(os.path.join("tests","TEST1","TEST1_align_position_correspondence.tsv")))
    
	align_position_correspondence = pd.read_csv(os.path.join("tests","TEST1","TEST1_align_position_correspondence.tsv"), sep='\t', header=None)
	assert(align_position_correspondence.iloc[9,3] == "9 (T)")
	assert(align_position_correspondence.iloc[26,11] == "25.1")
    
	assert(os.path.exists(os.path.join("tests","TEST1","TEST1_align_profile.fasta")))
	assert(os.path.exists(os.path.join("tests","TEST1","TEST1_align_profile.tsv")))
    
	align_profile = pd.read_csv(os.path.join("tests","TEST1","TEST1_align_profile.tsv"), sep='\t', header=None)
	assert(len(align_profile) == 15)
	assert(len(align_profile.columns) == 95)

	assert(os.path.exists(os.path.join("tests","TEST1","TEST1_clusterComposition.tsv")))
	cluster_composition = pd.read_csv(os.path.join("tests","TEST1","TEST1_clusterComposition.tsv"), sep='\t', header=None)
	assert(len(cluster_composition) == 207)
	assert(len(cluster_composition.columns) == 4)
    
	assert(os.path.exists(os.path.join("tests","TEST1","TEST1_country_origin_source_isolation_freq_matrix.tsv")))
	assert(os.path.exists(os.path.join("tests","TEST1","TEST1_country_origin_source_isolation_freq_pivot.tsv")))
	assert(os.path.exists(os.path.join("tests","TEST1","TEST1_dist.tsv")))
	assert(os.path.exists(os.path.join("tests","TEST1","TEST1_dist_grapetree.tsv")))
	assert(os.path.exists(os.path.join("tests","TEST1","TEST1_dist_hamming.tsv")))
    
	assert(os.path.exists(os.path.join("tests","TEST1","TEST1_metadata_w_partitions.tsv")))
	metadata_w_partitions = pd.read_csv(os.path.join("tests","TEST1","TEST1_metadata_w_partitions.tsv"), sep='\t', header=None)
	assert(metadata_w_partitions.iloc[0,9] == "nomenclature_code_" + "%s" % (datetime.date.today()))
	assert(metadata_w_partitions.iloc[0,11] == "MST-21x1.0") 
	assert(metadata_w_partitions.iloc[0,14] == "QUAL_called")
	assert(metadata_w_partitions.iloc[12,9] == "C2-C2-C2.1")        
    
	assert(os.path.exists(os.path.join("tests","TEST1","TEST1_metrics.tsv")))
    
	assert(os.path.exists(os.path.join("tests","TEST1","TEST1_nomenclature_changes.tsv")))
	nomenclature_changes = pd.read_csv(os.path.join("tests","TEST1","TEST1_nomenclature_changes.tsv"), sep='\t', header=None)
	assert(nomenclature_changes.iloc[1,2] == "7")
	assert(nomenclature_changes.iloc[9,6] == "1")

	assert(os.path.exists(os.path.join("tests","TEST1","TEST1_partitions.tsv")))    
	assert(os.path.exists(os.path.join("tests","TEST1","TEST1_partitions_summary.tsv")))
	partitions_summary = pd.read_csv(os.path.join("tests","TEST1","TEST1_partitions_summary.tsv"), sep='\t', header=None)
	assert(partitions_summary.iloc[1,9] ==  'A (83.3%), X Y Z (16.7%) (n = 12)') 

	assert(os.path.exists(os.path.join("tests","TEST1","TEST1_pos_report.tsv")))
	assert(os.path.exists(os.path.join("tests","TEST1","TEST1_redundantSamples.txt")))
	assert(os.path.exists(os.path.join("tests","TEST1","TEST1_sample9_closest5.zip")))
	assert(os.path.exists(os.path.join("tests","TEST1","TEST1_samples_excluded.txt")))
	assert(os.path.exists(os.path.join("tests","TEST1","TEST1_stableRegions.tsv")))
    
	shutil.rmtree(os.path.join("tests","TEST1"))
    
