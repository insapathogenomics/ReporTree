from partitioning_HC import HC
from datetime import datetime

print(datetime.now())
hc = HC('test_a')
hc.allele_profile = '/test_data/alleles.tsv'
hc.run()
print(datetime.now())