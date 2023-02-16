python /scripts/partitioning_HC.py -d_mx /test_data/distance_matrix.tsv -o d_single_direct
diff /tests/d_single/reference/test_d_clusterComposition.tsv $TMPDIR/d_single_direct/d_single_direct_clusterComposition.tsv
diff /tests/d_single/reference/test_d_partitions.tsv $TMPDIR/d_single_direct/d_single_direct_partitions.tsv
diff /tests/d_single/reference/test_d_single_HC.nwk $TMPDIR/d_single_direct/d_single_direct_single_HC.nwk