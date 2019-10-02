from cyvcf2 import VCF
from pyfaidx import Fasta
import sys
import os
import pandas as pd
import numpy as np
from variant import Variant
import itertools
from collections import Counter, defaultdict


def complement(c):
    base_pairs = {
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A',
        'N': 'N'
    }
    try:
        return base_pairs[c.upper()]
    except KeyError:
        raise ValueError(c + " is not a valid nucleotide.")


def get_complementary_sequence(sequence):
    """
    Returns a string of nucleotides complementary to the input string
    All letters in input sequence must A, C, T, or G, otherwise will raise a ValueError
    """
    comp_seq = []
    for c in sequence[::-1]:  # take the reverse complement
        comp_seq.append(complement(c))
    return "".join(comp_seq)


def generate_kmers(k):
    """Generates a list of all possible DNA sequences of length k. E.g. generate_kmers(2) will return
    [ AA, AC, AT, AG, CA, CC, CT, CG, TA, TC, TT, TG, GA, GC, GT, GG ] """
    len_k = int(k)
    if len_k < 0:
        raise ValueError("Must be a positive integer")
    combos = list(itertools.product('ACTG', repeat=len_k))
    seqs = []
    for seq in combos:
        seqs.append(''.join(seq))
    return seqs


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


def find_ref_kmer_freq(kmer_length):
    expected_path = 'data_files/chr22_' + str(kmer_length) + 'mer_frequency.csv'
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
            counts[next_seq] += 1
    print('Reference kmer frequency population done.')
    outfile = pd.DataFrame.from_dict(counts, orient='index')
    outfile.to_csv(expected_path)
    outfile = pd.read_csv(expected_path, header=None, skiprows=1, index_col=0, names=col_names)

    return outfile


def is_vcf(f_path):
    tokens = os.path.splitext(f_path)
    if tokens[1] == '.vcf' or tokens[1] == '.gz':
        return True
    return False


def is_quality_variant(var_to_test):
    """
    high quality variants will have FILTER == None
    Additionally, variants shoud be singletons ('AC' == 1) meaning that it is a unique observation
    AND we are ignoring insertions and deltions here
    """
    return var_to_test.FILTER is None and var_to_test.INFO.get('AC') == 1 and len(var_to_test.ALT) == 1 \
           and len(var_to_test.REF) == 1 and len(var_to_test.ALT[0]) == 1


def import_variants(f_path):
    return pd.read_csv(f_path, sep='\t')


def complete_sequence(adj_seq):
    """If 'N' is present iin the sequence, kmer is undefined"""
    return not ('N' in adj_seq or 'n' in adj_seq)


def process_variants(variants, kmer_size=7):
    """Input must be a pandas dataframe"""
    expected_path = 'data_files/variant_freq_chr22_7mers.csv'
    # if os.path.exists(expected_path):
    #     df = pd.read_csv(expected_path, header=None)
    #     return df
    if not kmer_size % 2 == 1:
        kmer_size += 1
    kmer_freq = find_ref_kmer_freq(kmer_size)  # frequency of 7mers on chr22
    fa = Fasta('/Users/simonelongo/too_big_for_icloud/REFERENCE_GENOME_GRch37.fa')
    chr22_ref = fa['22']
    mismatched_references = pd.DataFrame(columns=['POS', 'VCF_ref', 'FASTA_ref'])
    mismatch_count = 0
    counts = Counter()
    ref_nucleotides = defaultdict(list)
    start_idx_offset = int(kmer_size / 2 + 1)
    kmer_mid_idx = int(start_idx_offset - 1)  # also halfway index for kmer
    for i, r in variants.iterrows():
        position = r['POS']
        # take 7mer around variant. pyfaidx excludes start index and includes end index
        adj_seq = chr22_ref[(position - start_idx_offset):(position + kmer_mid_idx)].seq
        if complete_sequence(adj_seq):  # check if unknown nucleotides are present
            if adj_seq[kmer_mid_idx] != r['REF']:
                mismatched_references.iloc[mismatch_count] = [r['POS'], r['REF'], adj_seq[3]]
                mismatch_count += 1
            # now replace middle character with variant and append kmer_freq with associated sequence
            var_seq = list(adj_seq)
            var_seq[kmer_mid_idx] = r['ALT']
            new_seq = "".join(var_seq)
            counts[new_seq] += 1
            ref_nucleotides[new_seq].append(r['REF'])

    # convert lists of nucleotides to strings to handle more easily
    ref_nuc = defaultdict(str)
    for k, v in ref_nucleotides.items():
        ref_nuc[k] = "".join(v)

    var_count = pd.DataFrame.from_dict(counts, orient='index')
    var_count_path = "data_files/" + str(kmer_size) + "kmer_variant_counts.csv"
    var_count.to_csv(var_count_path)
    var_count = pd.read_csv(var_count_path, header=None, skiprows=1, index_col=0, names=['var_count'])

    ref_bp = pd.DataFrame.from_dict(ref_nuc, orient='index')
    ref_bp.to_csv("data_files/var_kmer_reference_bp.csv")
    ref_bp = pd.read_csv("data_files/var_kmer_reference_bp.csv", header=None, skiprows=1, index_col=0, names=['REF'])

    # join data
    joined_freq = kmer_freq.join(var_count, how='outer').join(ref_bp, how='outer')
    # save the files
    mismatched_references.to_csv('data_files/reference_mismatches.csv')
    joined_freq.to_csv("final_file.csv")
    return joined_freq


def impose_distance_requirement(sorted_vars, dist_bw_variants):
    # TODO
    print('No distance requirment imposed. Onwards!')
    return


if __name__ == "__main__":
    filename = sys.argv[1]
    variant_positions = defaultdict(Variant)
    kmer_len = 3
    if is_vcf(filename):  # import vcf file
        for variant in VCF(filename):
            if is_quality_variant(variant):
                # join is required because 'ALT' is returned as a list
                variant_positions[variant.POS] = Variant(variant.REF, "".join(variant.ALT), variant.POS)
        saved_csv_name = "chr22_variant_singletons.csv"
        generate_csv_from_variants(variant_positions, outfile=saved_csv_name)
        variant_singletons = import_variants(saved_csv_name)
    else:
        variant_singletons = import_variants(filename)
    print("Variants imported and saved.")
    if len(sys.argv) > 1:  # impose user supplied minimum distance between variants
        try:
            kmer_len = int(sys.argv[2])
            # dist_bw_variants = int(sys.argv[2])
            # if dist_bw_variants < 1:
            #     raise ValueError("Minimum distance must be positive integer!")
            # sorted_vars = variant_singletons.sort_values(by=['POS'])
            # sorted_vars['diff'] = sorted_vars['POS'].diff()
            # impose_distance_requirement(sorted_vars, dist_bw_variants)
        except ValueError:
            print("Invalid distance supplied. Running without minimum distance")
    process_variants(variant_singletons, kmer_size=kmer_len)
    # print(variant_singletons)
