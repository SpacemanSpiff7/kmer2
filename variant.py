class Variant:
    def __init__(self, ref, alt, pos):
        self.REF = ref
        self.ALT = alt
        self.POS = pos

    def __str__(self):
        # return "POS: " + str(self.POS) + "\tREF: " + str(self.REF) + "\tALT: " + str(self.ALT) + '\n'
        return str(self.POS) + "\t" + str(self.REF) + "\t" + str(self.ALT) + '\n'

    def __repr__(self):
        return "POS: " + str(self.POS) + "\tREF: " + str(self.REF) + "\tALT: " + str(self.ALT) + '\n'

    def csv_str(self):
        return str(self.POS) + "\t" + str(self.REF) + "\t" + str(self.ALT) + '\n'
