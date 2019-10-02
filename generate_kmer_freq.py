from pyfaidx import Fasta
from collections import Counter
import pandas as pd

kmer_length = 7
expected_path = 'data_files/chr22_7mer_frequency.csv'
counts = Counter()

ref_genome = Fasta('/Users/simonelongo/too_big_for_icloud/REFERENCE_GENOME_GRch37.fa')
ref_seq = str(ref_genome["22"])
for i in range(len(ref_seq) - (kmer_length - 1)):  # This takes the (1-based) reference sequence for chromosome 22
    next_seq = ref_seq[i:(i + kmer_length)]
    if not ('N' in next_seq or 'n' in next_seq):
        counts[next_seq] += 1
print('Reference kmer frequency population done.')

outfile = pd.DataFrame.from_dict(counts, orient='index')
outfile.to_csv(expected_path)
