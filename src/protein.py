from Bio.Seq import Seq


def translate_rna(rna):
    return str(Seq(rna).translate())

def codon_usage(rna):
    codons = {}

    for i in range(0, len(rna) - 2, 3):
        codon = rna[i:i+3]

        if codon in codons:
            codons[codon] += 1
        else:
            codons[codon] = 1

    return codons

def amino_acid_composition(protein):
    composition = {}

    for amino_acid in protein:
        if amino_acid in composition:
            composition[amino_acid] += 1
        else:
            composition[amino_acid] = 1

    return composition

AMINO_ACID_MASS = {
    "A": 89.09,
    "R": 174.20,
    "N": 132.12,
    "D": 133.10,
    "C": 121.15,
    "Q": 146.15,
    "E": 147.13,
    "G": 75.07,
    "H": 155.16,
    "I": 131.17,
    "L": 131.17,
    "K": 146.19,
    "M": 149.21,
    "F": 165.19,
    "P": 115.13,
    "S": 105.09,
    "T": 119.12,
    "W": 204.23,
    "Y": 181.19,
    "V": 117.15,
}

def molecular_weight(protein):
    weight = 0

    for aa in protein:
        if aa in AMINO_ACID_MASS:
            weight += AMINO_ACID_MASS[aa]

    return round(weight, 2)

def codon_usage(sequence):
    usage = {}

    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]

        if len(codon) == 3:
            usage[codon] = usage.get(codon, 0) + 1

    return usage