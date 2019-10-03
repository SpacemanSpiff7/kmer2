from collections import defaultdict
from cyvcf2 import VCF
from kmer import is_quality_variant
from variant import Variant

filename = '../samp.vcf'
variant_positions = defaultdict(Variant)
for variant in VCF(filename):
    if is_quality_variant(variant):
        # join is required because 'ALT' is returned as a list
        variant_positions[variant.POS] = Variant(variant.REF, "".join(variant.ALT), variant.POS)


print(variant_positions)
