from datetime import datetime
import sys
sys.path.append('/scripts/')

from partitioning_HC import HC

print(datetime.now())
hc = HC('a_single_as_library')
hc.out = 'a_single_as_library'
hc.allele_profile = '/test_data/alleles.tsv'
hc.run()
print(datetime.now())