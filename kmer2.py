import sys
from collections import defaultdict, Counter
from cyvcf2 import VCF
from pyfaidx import Fasta
from kmer import is_vcf, find_ref_kmer_freq, complete_sequence
from kmer_shift import is_quality_variant, generate_csv_from_variants
from variant import Variant
import pandas as pd


def process_variants(variants, kmer_size):
    if not kmer_size % 2 == 1:
        kmer_size += 1  # kmer size must be an odd number
    kmer_freq = find_ref_kmer_freq(kmer_size)
    fa = Fasta('/Users/simonelongo/too_big_for_icloud/REFERENCE_GENOME_GRch37.fa')
    chr22_ref = fa['22']  # reference sequence of chr22
    transitions = defaultdict(Counter)
    start_idx_offset = int(kmer_size / 2 + 1)
    kmer_mid_idx = int(start_idx_offset - 1)  # also halfway index for kmer
    for idx, row in variants.iterrows():
        position = row['POS']
        # take 7mer around variant. pyfaidx excludes start index and includes end index
        adj_seq = chr22_ref[(position - start_idx_offset):(position + kmer_mid_idx)].seq
        if complete_sequence(adj_seq):
            transitions[adj_seq][row['ALT']] += 1
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
    process_variants(variant_singletons, kmer_size=kmer_len)
