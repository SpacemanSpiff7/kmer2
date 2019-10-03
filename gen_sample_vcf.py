from os import system
import random
from pyfaidx import Fasta
from cyvcf2 import VCF, Writer


def gen_com(chrom, start, stop):
    return "tabix /Users/simonelongo/too_big_for_icloud/gnomad.genomes.r2.1.1.sites.vcf.bgz " + str(chrom) + ":" + str(start) + "-" + str(
        stop) + " >> samp.vcf"


fa = Fasta('/Users/simonelongo/too_big_for_icloud/REFERENCE_GENOME_GRch37.fa')
vcf = VCF('/Users/simonelongo/too_big_for_icloud/gnomad.genomes.r2.1.1.sites.vcf.bgz')
write = Writer('samp.vcf', vcf)
write.write_header()
keys = list(fa.keys())

for key in keys[0:25]:
    if key != 'MT':
        begin = random.randint(1000, len(fa[key]) - 1000)
        system(gen_com(key, begin, begin + 1000))

write.close()
