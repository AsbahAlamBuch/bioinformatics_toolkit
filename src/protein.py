from Bio.Seq import Seq


def translate_rna(rna):
    return str(Seq(rna).translate())