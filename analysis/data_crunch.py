"""
Created by Simone Longo Sept 24, 2019
"""
import pandas as pd
import numpy as np
import statistics
import matplotlib.pyplot as plt
from collections import Counter, defaultdict


# Final file used here does not have a label for the first column, which will be used as the index

def most_common_char(sequence, freq=False):
    split = Counter()
    for c in sequence:
        split[c] += 1

    max_nucleotide_count = 0
    most_freq_nuc = 'N'
    for c in split:
        if split[c] > max_nucleotide_count:
            max_nucleotide_count = split[c]
            most_freq_nuc = c
        # max_nucleotide_count = split[c] if split[c] > max_nucleotide_count else max_nucleotide_count
    if freq is True:
        return max_nucleotide_count / len(sequence)
    return most_freq_nuc


def common_char(sequence):
    return Counter(sequence).most_common(1)[0][0]


def char_freq(sequence):
    return str(Counter(sequence).most_common(1)[0][1] / len(sequence))


data = pd.read_csv('final_file_noidx.csv')
data['ALT'] = data['kmer'].str[3]
data['freq'] = data['var_count'] / data['ref_count']
data.sort_values(by=['freq'])

# take the most common nucleotide in the string and also calculate its frequency
data['most_common_ref'] = data['REF'].map(common_char, na_action='ignore')
data['most_common_ref_frac'] = data['REF'].map(char_freq, na_action='ignore')
"""
The above process takes a file containing the frequency of each 7mer on chr22 in addition to
Using variants from gnomad to find the frequency of each mutated 7mer.
It also looks at what the reference nucleotide is for each variant and stores that.
"""
data.to_csv("chr22_7mer_variants.csv")

# garbage visualization
plt.hist(data['freq'].to_list(), bins=50)
plt.show()

mutations = defaultdict(list)
for index, row in data.iterrows():
    if not pd.isnull(row['most_common_ref']):
        # store transition of tuple in format (reference nucleotide, mutated)
        transition = (row['most_common_ref'], row['ALT'])
        mutations[transition].append(row['freq'])  # append freq.
        """
        freq represents the occurrence of a 7mer with a mutation divided by the occurrence of that 7mer 
        in the reference sequence of chr22
        """

tr_fr = defaultdict(list)
for key, value in mutations.items():
    tr_fr[key] = [len(value), statistics.median(value)]

fr = pd.read_csv("frequency_mut_chr22.csv", skiprows=1, header=None)
fr.columns = ['transition', 'count', 'avg_freq']
fr.sort_values(by=['count'])
print(fr)
plt.bar(fr['transition'].to_list(), fr['count'].to_list())
plt.xticks(rotation=90)
plt.ylabel("Number of occurrences")
plt.xlabel("Transition")
plt.show()
#print(pd.DataFrame.from_dict(tr_fr, orient='index').to_csv("frequency_mut_chr22.csv"))
#print(tr_fr)
# print(data)
