python /scripts/partitioning_HC.py -a /test_data/alleles.tsv -o test_a
diff /tests/a_single/reference/test_a_clusterComposition.tsv test_a_clusterComposition.tsv
diff /tests/a_single/reference/test_a_dist.tsv test_a_dist.tsv
diff /tests/a_single/reference/test_a_partitions.tsv test_a_partitions.tsv
diff /tests/a_single/reference/test_a_single_HC.nwk test_a_single_HC.nwk