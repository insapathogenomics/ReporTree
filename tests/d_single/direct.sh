python /scripts/partitioning_HC.py -d_mx /test_data/distance_matrix.tsv -o test_d
diff /tests/d_single/reference/test_d_clusterComposition.tsv test_d_clusterComposition.tsv
diff /tests/d_single/reference/test_d_partitions.tsv test_d_partitions.tsv
diff /tests/d_single/reference/test_d_single_HC.nwk test_d_single_HC.nwk