from partitioning_HC import HC

hc = HC('test_a')
hc.allele_profile = '/test_data/alleles.tsv'
print(f"Output folder: {hc.out}")
print(f"Distance matrix: {hc.distance_matrix}")
print(f"Allele profile: {hc.allele_profile}")
print(f"Method threshold: {hc.method_threshold}")
print(f"Percentage method threshold: {hc.pct_HCmethod_threshold}")
print(f"Samples called: {hc.samples_called}")
print(f"Loci called: {hc.loci_called}")
print(f"Metadata: {hc.metadata}")
print(f"Filter column: {hc.filter_column}")
print(f"Distances: {hc.dist}")
hc.run()