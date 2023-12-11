"""
Runs Several Tests for ReporTree functionality
author daniel.sobral@insa.min-saude.pt
"""
import os
import shutil
import subprocess
import pandas as pd
import datetime

def test_HC():
	
	if not os.path.exists(os.path.join("tests","TEST3")):
		os.mkdir(os.path.join("tests","TEST3"))

	result = subprocess.check_output("python reportree.py " + \
		"-a tests/allele_matrix.tsv -out tests/TEST3/TEST3 --analysis HC -m tests/metadata.tsv " + \
		"-f \'country_origin != B;date >= 2023-02-02\' --subset --missing-code 0 --loci-called 0.98 " + \
        "--sample-ATCG-content  0.99 --nomenclature-file tests/nomenclature.tsv --nomenclature-code-levels single-8,country_origin " + \
        "--columns_summary_report first_seq_date,last_seq_date,timespan_days,country_origin,n_country_origin,source_isolation,n_source_isolation " + \
        "--partitions2report single-1-10", shell=True)
       

	assert(os.path.exists(os.path.join("tests","TEST3","TEST3_clusterComposition.tsv")))  
	assert(os.path.exists(os.path.join("tests","TEST3","TEST3_dist_hamming.tsv")))
	assert(os.path.exists(os.path.join("tests","TEST3","TEST3_flt_samples_matrix.tsv")))
	assert(os.path.exists(os.path.join("tests","TEST3","TEST3_loci_report.tsv")))
	assert(os.path.exists(os.path.join("tests","TEST3","TEST3_metadata_w_partitions.tsv"))) 
	metadata_w_partitions = pd.read_csv(os.path.join("tests","TEST3","TEST3_metadata_w_partitions.tsv"), sep='\t', header=None)
	assert(metadata_w_partitions.iloc[3,9] == "C1-X_Y_Z") 
	assert(os.path.exists(os.path.join("tests","TEST3","TEST3_nomenclature_changes.tsv")))
	current_date = "%s" % (datetime.date.today())
	assert(os.path.exists(os.path.join("tests","TEST3","TEST3_nomenclature_code_"+current_date+"_summary.tsv")))
	assert(os.path.exists(os.path.join("tests","TEST3","TEST3_partitions.tsv")))
	assert(os.path.exists(os.path.join("tests","TEST3","TEST3_partitions_summary.tsv")))
	assert(os.path.exists(os.path.join("tests","TEST3","TEST3_single_HC.nwk")))
	assert(os.path.exists(os.path.join("tests","TEST3","TEST3_subset_matrix.tsv")))
    
	shutil.rmtree(os.path.join("tests","TEST3"))
    