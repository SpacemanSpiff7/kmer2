class Variant:
    def __init__(self, ref, alt, pos, chrom=22):
        self.REF = ref
        self.ALT = alt
        self.POS = pos
        self.chrom = chrom

    def __str__(self):
        # return "POS: " + str(self.POS) + "\tREF: " + str(self.REF) + "\tALT: " + str(self.ALT) + '\n'
        return str(self.chrom) + "\t" + str(self.POS) + "\t" + str(self.REF) + "\t" + str(self.ALT) + '\n'

    def __repr__(self):
        return "CHROM: " + str(self.chrom) + "\t" + "POS: " + str(self.POS) + "\tREF: " + str(self.REF) + "\tALT: " + str(self.ALT) + '\n'

    def csv_str(self):
        return str(self.chrom) + "\t" + str(self.POS) + "\t" + str(self.REF) + "\t" + str(self.ALT) + '\n'
