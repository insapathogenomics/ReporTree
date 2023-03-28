"""
Runs Several Tests for ReporTree functionality
author daniel.sobral@insa.min-saude.pt
"""
import os
import shutil
import subprocess
import pandas as pd
import datetime

def test_treecluster():
	
	if not os.path.exists(os.path.join("tests","TEST4")):
		os.mkdir(os.path.join("tests","TEST4"))

	result = subprocess.check_output("python reportree.py -t tests/tree.nwk " + \
		"-out tests/TEST4/TEST4 --analysis HC -m tests/metadata.tsv --nomenclature-file tests/nomenclature.tsv " + \
		"--columns_summary_report first_seq_date,last_seq_date,timespan_days,country_origin,n_country_origin,source_isolation,n_source_isolation", shell=True)
      
	assert(os.path.exists(os.path.join("tests","TEST4","TEST4.log")))     
	assert(os.path.exists(os.path.join("tests","TEST4","TEST4_clusterComposition.tsv")))
	assert(os.path.exists(os.path.join("tests","TEST4","TEST4_metadata_w_partitions.tsv")))
    
	assert(os.path.exists(os.path.join("tests","TEST4","TEST4_nomenclature_changes.tsv")))
	nomenclature_changes = pd.read_csv(os.path.join("tests","TEST4","TEST4_nomenclature_changes.tsv"), sep='\t', header=None)
	assert(nomenclature_changes.iloc[2,2] == "multiple clusters")
	assert(nomenclature_changes.iloc[2,3] == "new (split_merge)")   
    
	assert(os.path.exists(os.path.join("tests","TEST4","TEST4_partitions.tsv")))
	assert(os.path.exists(os.path.join("tests","TEST4","TEST4_partitions_summary.tsv")))
        
	shutil.rmtree(os.path.join("tests","TEST4"))
    
