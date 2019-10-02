import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict

transitions = pd.read_csv('bp_counts_per3mer.csv', index_col=0)
kmer_freq = pd.read_csv('data_files/chr22_3mer_frequency.csv', index_col=0)
kmer_freq.columns = ['count']
merged = transitions.join(kmer_freq, how='inner')

kmers_idx = []
for seq in merged.index:
    if seq[1] == 'A' or seq[1] == 'T':
        kmers_idx.append(seq)

trunc = merged[merged.index.isin(kmers_idx)]
rel_trunc = trunc.iloc[:, 0:4].div(trunc['count'], axis=0)

flank_count = pd.DataFrame()
for idx, row in rel_trunc.iterrows():
    if idx[1] == 'A':
        flank_count.loc[idx[0], idx[2]] = row['G']
    if idx[1] == 'T':
        flank_count.loc[idx[0], idx[2]] = row['C']

flank_count = flank_count.sort_index(ascending=False)
flank_count = flank_count.reindex(sorted(flank_count.columns), axis=1)
plt.matshow(flank_count)
plt.colorbar()
plt.title('A > G Transitions of 3-mers')
plt.xticks(range(4), list(flank_count.columns))
plt.yticks(range(4), list(flank_count.index))
# relative_freq = merged.iloc[:, 0:6].div(merged['count'], axis=0)
#
# flank_count = pd.DataFrame()
# # considering A > G transitions to compare to Carlson paper
# for idx, row in relative_freq.iterrows():
#     if pd.notna(row['AG']):
#         flank_count.loc[idx[0], idx[2]] = row['AG']
# print(flank_count)
# plt.matshow(flank_count)

# plt.show()
