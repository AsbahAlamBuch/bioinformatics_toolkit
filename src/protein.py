from Bio.Seq import Seq


def translate_rna(rna):
    return str(Seq(rna).translate())

def codon_usage(sequence):
    codons = {}

    for i in range(0, len(sequence) - 2, 3):

        codon = sequence[i:i + 3]

        if len(codon) == 3:

            if codon in codons:
                codons[codon] += 1
            else:
                codons[codon] = 1

    return codons