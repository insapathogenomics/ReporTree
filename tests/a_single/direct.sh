python /scripts/partitioning_HC.py -a /test_data/alleles.tsv -o a_single_direct
diff /tests/a_single/reference/test_a_clusterComposition.tsv $TMPDIR/a_single_direct/a_single_direct_clusterComposition.tsv
diff /tests/a_single/reference/test_a_dist.tsv $TMPDIR/a_single_direct/a_single_direct_dist.tsv
diff /tests/a_single/reference/test_a_partitions.tsv  $TMPDIR/a_single_direct/a_single_direct_partitions.tsv
diff /tests/a_single/reference/test_a_single_HC.nwk  $TMPDIR/a_single_direct/a_single_direct_single_HC.nwk
