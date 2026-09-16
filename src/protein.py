from typing import Dict

from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis


VALID_AMINO_ACIDS = set(
    "ARNDCQEGHILKMFPSTWYV"
)


def validate_protein(protein: str) -> bool:
    """Check whether a protein sequence contains valid amino acids."""

    if not protein:
        return False

    protein = protein.upper()

    return all(
        amino_acid in VALID_AMINO_ACIDS
        for amino_acid in protein
    )


def translate_rna(rna: str) -> str:
    """Translate an RNA sequence into a protein sequence."""

    if not rna:
        return ""

    rna = rna.upper()

    if any(base not in "AUGC" for base in rna):
        raise ValueError("Invalid RNA sequence.")

    # Ignore incomplete trailing nucleotides.
    complete_length = len(rna) - (len(rna) % 3)

    rna = rna[:complete_length]

    if not rna:
        return ""

    return str(
        Seq(rna).translate(
            to_stop=False
        )
    )


def amino_acid_composition(
    protein: str
) -> Dict[str, int]:
    """
    Calculate amino acid composition.

    Only amino acids actually present in the sequence
    are included in the returned dictionary.
    """

    if not validate_protein(protein):
        raise ValueError("Invalid protein sequence.")

    protein = protein.upper()

    composition = {}

    for amino_acid in protein:

        composition[amino_acid] = (
            composition.get(amino_acid, 0) + 1
        )

    return composition


def molecular_weight(
    protein: str
) -> float:
    """Calculate molecular weight of a protein in Daltons."""

    if not validate_protein(protein):
        raise ValueError("Invalid protein sequence.")

    analysis = ProteinAnalysis(
        protein.upper()
    )

    return round(
        analysis.molecular_weight(),
        2
    )


def calculate_gravy(
    protein: str
) -> float:
    """Calculate the GRAVY hydropathy score."""

    if not validate_protein(protein):
        raise ValueError("Invalid protein sequence.")

    analysis = ProteinAnalysis(
        protein.upper()
    )

    return round(
        analysis.gravy(),
        3
    )


def calculate_pI(
    protein: str
) -> float:
    """Calculate the theoretical isoelectric point."""

    if not validate_protein(protein):
        raise ValueError("Invalid protein sequence.")

    analysis = ProteinAnalysis(
        protein.upper()
    )

    return round(
        analysis.isoelectric_point(),
        2
    )


def get_protein_length(
    protein: str
) -> int:
    """Return the length of a valid protein sequence."""

    if not validate_protein(protein):
        raise ValueError("Invalid protein sequence.")

    return len(protein)


def calculate_instability_index(
    protein: str
) -> float:
    """Calculate the protein instability index."""

    if not validate_protein(protein):
        raise ValueError("Invalid protein sequence.")

    analysis = ProteinAnalysis(
        protein.upper()
    )

    return round(
        analysis.instability_index(),
        2
    )


def codon_usage(
    rna: str
) -> Dict[str, int]:
    """
    Count codon usage in an RNA sequence.

    Incomplete trailing nucleotides are ignored.
    """

    if not rna:
        return {}

    rna = rna.upper()

    if any(base not in "AUGC" for base in rna):
        raise ValueError("Invalid RNA sequence.")

    usage = {}

    for i in range(
        0,
        len(rna) - 2,
        3
    ):

        codon = rna[i:i + 3]

        usage[codon] = usage.get(
            codon,
            0
        ) + 1

    return usage