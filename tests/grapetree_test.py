"""
Runs Several Tests for ReporTree functionality
author daniel.sobral@insa.min-saude.pt
"""
import os
import shutil
import subprocess
import pandas as pd
import datetime

def test_grapetree():
	
	if not os.path.exists(os.path.join("tests","TEST2")):
		os.mkdir(os.path.join("tests","TEST2"))

	result = subprocess.check_output("python reportree.py " + \
		"-a tests/allele_matrix.tsv -out tests/TEST2/TEST2 --analysis grapetree -m tests/metadata.tsv " + \
		"-f \'country_origin != B;date >= 2023-02-02\' --subset --missing-code 0 --loci-called 0.98 " + \
        "--sample-ATCG-content  0.99 --nomenclature-file tests/nomenclature.tsv --nomenclature-code-levels 8,country_origin " + \
        "--columns_summary_report first_seq_date,last_seq_date,timespan_days,country_origin,n_country_origin,source_isolation,n_source_isolation " + \
        "--partitions2report 1-10", shell=True)
       

    ### Check if files exist:
	assert(os.path.exists(os.path.join("tests","TEST2","TEST2.log")))
	assert(os.path.exists(os.path.join("tests","TEST2","TEST2.nwk")))
	assert(os.path.exists(os.path.join("tests","TEST2","TEST2_clusterComposition.tsv")))
	assert(os.path.exists(os.path.join("tests","TEST2","TEST2_dist.tsv")))
	assert(os.path.exists(os.path.join("tests","TEST2","TEST2_dist_grapetree.tsv")))
	assert(os.path.exists(os.path.join("tests","TEST2","TEST2_dist_hamming.tsv")))
	assert(os.path.exists(os.path.join("tests","TEST2","TEST2_flt_samples_matrix.tsv")))
	assert(os.path.exists(os.path.join("tests","TEST2","TEST2_loci_report.tsv")))

	assert(os.path.exists(os.path.join("tests","TEST2","TEST2_metadata_w_partitions.tsv")))
	metadata_w_partitions = pd.read_csv(os.path.join("tests","TEST2","TEST2_metadata_w_partitions.tsv"), sep='\t', header=None)
	assert(metadata_w_partitions.iloc[3,9] == "C1-X_Y_Z") 

	assert(os.path.exists(os.path.join("tests","TEST2","TEST2_nomenclature_changes.tsv")))
	current_date = "%s" % (datetime.date.today())
	assert(os.path.exists(os.path.join("tests","TEST2","TEST2_nomenclature_code_"+current_date+"_summary.tsv")))
	assert(os.path.exists(os.path.join("tests","TEST2","TEST2_partitions.tsv")))
	assert(os.path.exists(os.path.join("tests","TEST2","TEST2_partitions_summary.tsv")))
	assert(os.path.exists(os.path.join("tests","TEST2","TEST2_redundantSamples.txt")))
	assert(os.path.exists(os.path.join("tests","TEST2","TEST2_subset_matrix.tsv")))

	shutil.rmtree(os.path.join("tests","TEST2"))
    