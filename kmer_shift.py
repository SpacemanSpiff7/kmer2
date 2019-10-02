from cyvcf2 import VCF
from pyfaidx import Fasta
import sys
import os
import pandas as pd
import numpy as np

from kmer import complete_sequence, complement
from kmer_seq import Kmer
from variant import Variant
import itertools
from collections import Counter, defaultdict


def find_ref_kmer_freq2(kmer_length):
    expected_path = 'data_files/chr22_' + str(kmer_length) + 'mer_frequency_comp.csv'
    col_names = ['ref_count']
    if os.path.exists(expected_path):
        df = pd.read_csv(expected_path, header=None, index_col=0, names=col_names)
        print('Reference kmer frequency successfully read.')
        return df
    counts = Counter()
    ref_genome = Fasta('/Users/simonelongo/too_big_for_icloud/REFERENCE_GENOME_GRch37.fa')
    ref_seq = str(ref_genome["22"])
    for i in range(len(ref_seq) - (kmer_length - 1)):  # This takes the (1-based) reference sequence for chromosome 22
        next_seq = ref_seq[i:(i + kmer_length)]
        if not ('N' in next_seq or 'n' in next_seq):
            counts[Kmer(next_seq)] += 1
    print('Reference kmer frequency population done.')
    outfile = pd.DataFrame.from_dict(counts, orient='index')
    outfile.to_csv(expected_path)
    outfile = pd.read_csv(expected_path, header=None, skiprows=1, index_col=0, names=col_names)

    return outfile


def is_quality_variant(var_to_test):
    """
    high quality variants will have FILTER == None
    Additionally, variants shoud be singletons ('AC' == 1) meaning that it is a unique observation
    AND we are ignoring insertions and deltions here
    """
    return var_to_test.FILTER is None and var_to_test.INFO.get('AC') == 1 and len(var_to_test.ALT) == 1 \
           and len(var_to_test.REF) == 1 and len(var_to_test.ALT[0]) == 1


def is_vcf(f_path):
    tokens = os.path.splitext(f_path)
    if tokens[1] == '.vcf' or tokens[1] == '.gz':
        return True
    return False


def generate_csv_from_variants(variants, outfile="variants.csv"):
    """Converts a dictionary to a csv to prevemt redundant slow operations"""
    if not type(variants) == dict:
        print("Input must be a dictionary")
        return
    output = open(outfile, "w+")
    output.write("POS\tREF\tALT\n")  # header same for all
    for k, v in variants.items():
        output.write(str(v))
    output.close()


def get_transition(ref, alt):
    if ref != 'A' and ref != 'C':
        ref = complement(ref)
        alt = complement(alt)
    return ref + alt


def process_variants2(variants, kmer_size):
    if not kmer_size % 2 == 1:
        kmer_size += 1  # kmer size must be an odd number
    fa = Fasta('/Users/simonelongo/too_big_for_icloud/REFERENCE_GENOME_GRch37.fa')
    chr22_ref = fa['22']  # reference sequence of chr22
    # TODO find ref kmer freq considering complement
    kmer_freq = find_ref_kmer_freq2(kmer_size)
    transitions = defaultdict(Counter)
    start_idx_offset = int(kmer_size / 2 + 1)
    kmer_mid_idx = int(start_idx_offset - 1)  # also halfway index for kmer
    for idx, row in variants.iterrows():
        position = row['POS']
        # take 7mer around variant. pyfaidx excludes start index and includes end index
        adj_seq = chr22_ref[(position - start_idx_offset):(position + kmer_mid_idx)].seq
        if complete_sequence(adj_seq):
            transitions[Kmer(adj_seq)][get_transition(adj_seq[kmer_mid_idx], row['ALT'])] += 1
    count_fpath = "bp_counts_per" + str(kmer_len) + "mer.csv"
    transitions_df = pd.DataFrame.from_dict(transitions, orient='index')
    merged_df = pd.merge(transitions_df, kmer_freq, left_index=True, right_index=True, how='inner')
    transitions_df.to_csv(count_fpath)


def process_variants(variants, kmer_size):
    if not kmer_size % 2 == 1:
        kmer_size += 1  # kmer size must be an odd number
    from kmer import find_ref_kmer_freq
    kmer_freq = find_ref_kmer_freq(kmer_size)
    fa = Fasta('/Users/simonelongo/too_big_for_icloud/REFERENCE_GENOME_GRch37.fa')
    chr22_ref = fa['22']  # reference sequence of chr22
    mismatched_references = pd.DataFrame(columns=['POS', 'VCF_ref', 'FASTA_ref'])
    mismatch_count = 0
    counts = Counter()
    transitions = defaultdict(Counter)
    ref_nucleotides = defaultdict(list)
    # indices for genome positions
    start_idx_offset = int(kmer_size / 2 + 1)
    kmer_mid_idx = int(start_idx_offset - 1)  # also halfway index for kmer
    for idx, row in variants.iterrows():
        position = row['POS']
        # take 7mer around variant. pyfaidx excludes start index and includes end index
        adj_seq = chr22_ref[(position - start_idx_offset):(position + kmer_mid_idx)].seq
        if complete_sequence(adj_seq):
            if adj_seq[kmer_mid_idx] != row['REF']:
                mismatched_references.iloc[mismatch_count] = [row['POS'], row['REF'], adj_seq[3]]
                mismatch_count += 1
            # Replace middle ref nucleotide with the ALT
            var_seq = list(adj_seq)
            var_seq[kmer_mid_idx] = row['ALT']
            new_seq = "".join(var_seq)
            counts[new_seq] += 1
            transitions[Kmer(adj_seq)][row['ALT']] += 1
            ref_nucleotides[new_seq].append(row['REF'])
    # convert lists of nucleotides to strings to handle more easily
    # ref_nuc = defaultdict(str)
    # for k, v in ref_nucleotides.items():
    #     ref_nuc[k] = "".join(v)
    count_fpath = "bp_counts_per" + str(kmer_len) + "mer.csv"
    transitions_df = pd.DataFrame.from_dict(transitions, orient='index')
    transitions_df.to_csv(count_fpath)


if __name__ == "__main__":
    filename = sys.argv[1]
    variant_positions = defaultdict(Variant)
    kmer_len = 3  # this is the default
    if is_vcf(filename):  # import vcf file
        for variant in VCF(filename):
            if is_quality_variant(variant):
                # join is required because 'ALT' is returned as a list
                variant_positions[variant.POS] = Variant(variant.REF, "".join(variant.ALT), variant.POS)
        saved_csv_name = "chr22_variant_singletons.csv"
        generate_csv_from_variants(variant_positions, outfile=saved_csv_name)
        variant_singletons = pd.read_csv(saved_csv_name, sep='\t')
    else:
        variant_singletons = pd.read_csv(filename, sep='\t')
    print("Variants imported and saved.")
    if len(sys.argv) > 2:  # kmer size from user input
        try:
            kmer_len = int(sys.argv[2])
        except:
            print("Invalid distance supplied. Running with default kmer size of " + str(kmer_len))
    process_variants2(variant_singletons, kmer_size=kmer_len)
