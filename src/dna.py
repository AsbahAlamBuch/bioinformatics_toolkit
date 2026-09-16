from typing import List, Tuple, Dict
from Bio.SeqUtils.MeltingTemp import Tm_NN


def get_sequence_length(sequence: str) -> int:
    """Return the length of a DNA sequence."""
    return len(sequence)


def validate_dna(sequence: str) -> bool:
    """Check whether a sequence contains only valid DNA bases."""

    if not sequence:
        return False

    sequence = sequence.upper()

    return all(
        base in "ATGC"
        for base in sequence
    )


def calculate_gc_content(sequence: str) -> float:
    """Calculate the GC content of a valid DNA sequence as a percentage."""

    if not validate_dna(sequence):
        raise ValueError("Invalid DNA sequence.")

    sequence = sequence.upper()

    gc_count = sequence.count("G") + sequence.count("C")

    return (gc_count / len(sequence)) * 100


def count_nucleotides(sequence: str) -> Dict[str, int]:
    """Count the number of A, T, G and C nucleotides."""

    if not validate_dna(sequence):
        raise ValueError("Invalid DNA sequence.")

    sequence = sequence.upper()

    counts = {
        "A": sequence.count("A"),
        "T": sequence.count("T"),
        "G": sequence.count("G"),
        "C": sequence.count("C"),
    }

    return counts


def reverse_sequence(sequence: str) -> str:
    """Return the reverse of a valid DNA sequence."""

    if not validate_dna(sequence):
        raise ValueError("Invalid DNA sequence.")

    return sequence.upper()[::-1]


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a DNA sequence."""

    if not validate_dna(sequence):
        raise ValueError("Invalid DNA sequence.")

    complement = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G",
    }

    sequence = sequence.upper()

    return "".join(
        complement[base]
        for base in reversed(sequence)
    )


def transcribe_dna(sequence: str) -> str:
    """Transcribe a valid DNA sequence into RNA."""

    if not validate_dna(sequence):
        raise ValueError("Invalid DNA sequence.")

    sequence = sequence.upper()

    return sequence.replace("T", "U")


def find_motif(sequence: str, motif: str) -> List[int]:
    """
    Find all occurrences of a motif in a DNA sequence.

    Positions are returned using 0-based indexing.
    """

    if not validate_dna(sequence):
        raise ValueError("Invalid DNA sequence.")

    if not motif:
        raise ValueError("Motif cannot be empty.")

    motif = motif.upper()

    if not validate_dna(motif):
        raise ValueError("Invalid DNA motif.")

    sequence = sequence.upper()

    positions = []

    for i in range(len(sequence) - len(motif) + 1):

        if sequence[i:i + len(motif)] == motif:
            positions.append(i)

    return positions


def find_orfs(
    sequence: str,
) -> List[Tuple[int, int, str]]:
    """
    Find open reading frames in all six reading frames.

    Searches:
    - Three reading frames on the original strand.
    - Three reading frames on the reverse-complement strand.

    An ORF begins with ATG and ends at the first in-frame
    stop codon (TAA, TAG, or TGA).

    Returns
    -------
    list of tuples
        Each tuple contains:
        (start_position, end_position, ORF_sequence)

    Positions are mapped to the original DNA sequence
    coordinates and use 0-based indexing.

    For reverse-strand ORFs, start_position will be greater
    than end_position because the ORF runs in the opposite
    direction.
    """

    if not validate_dna(sequence):
        raise ValueError("Invalid DNA sequence.")

    sequence = sequence.upper()

    orfs = []

    sequence_length = len(sequence)

    # Forward strand

    for frame in range(3):

        start = frame

        while start <= sequence_length - 3:

            if sequence[start:start + 3] == "ATG":

                for end in range(
                    start + 3,
                    sequence_length - 2,
                    3
                ):

                    codon = sequence[end:end + 3]

                    if codon in {
                        "TAA",
                        "TAG",
                        "TGA"
                    }:

                        orf_sequence = sequence[
                            start:end + 3
                        ]

                        orfs.append(
                            (
                                start,
                                end + 3,
                                orf_sequence
                            )
                        )

                        break

            start += 3

    # Reverse-complement strand

    reverse_sequence_data = reverse_complement(sequence)

    for frame in range(3):

        start = frame

        while start <= sequence_length - 3:

            if reverse_sequence_data[start:start + 3] == "ATG":

                for end in range(
                    start + 3,
                    sequence_length - 2,
                    3
                ):

                    codon = reverse_sequence_data[
                        end:end + 3
                    ]

                    if codon in {
                        "TAA",
                        "TAG",
                        "TGA"
                    }:

                        orf_sequence = reverse_sequence_data[
                            start:end + 3
                        ]

                        original_start = (
                            sequence_length - start
                        )

                        original_end = (
                            sequence_length - (end + 3)
                        )

                        orfs.append(
                            (
                                original_start,
                                original_end,
                                orf_sequence
                            )
                        )

                        break

            start += 3

    return orfs


def kmer_count(
    sequence: str,
    k: int
) -> Dict[str, int]:
    """Count all k-mers in a valid DNA sequence."""

    if not validate_dna(sequence):
        raise ValueError("Invalid DNA sequence.")

    if k <= 0:
        raise ValueError("k must be greater than 0.")

    if k > len(sequence):
        return {}

    sequence = sequence.upper()

    kmers = {}

    for i in range(len(sequence) - k + 1):

        kmer = sequence[i:i + k]

        kmers[kmer] = kmers.get(kmer, 0) + 1

    return kmers


RESTRICTION_ENZYMES = {
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
    "HindIII": "AAGCTT",
    "PstI": "CTGCAG",
    "SmaI": "CCCGGG",
    "XhoI": "CTCGAG",
}


def find_restriction_sites(
    sequence: str
) -> Dict[str, List[int]]:
    """
    Find common restriction enzyme recognition sites
    in a valid DNA sequence.

    Positions are returned using 0-based indexing.
    """

    if not validate_dna(sequence):
        raise ValueError("Invalid DNA sequence.")

    sequence = sequence.upper()

    results = {}

    for enzyme, recognition_site in RESTRICTION_ENZYMES.items():

        positions = []

        for i in range(
            len(sequence) - len(recognition_site) + 1
        ):

            if sequence[
                i:i + len(recognition_site)
            ] == recognition_site:

                positions.append(i)

        if positions:
            results[enzyme] = positions

    return results


def calculate_melting_temperature(sequence: str) -> float:
    """
    Calculate an approximate DNA melting temperature (Tm).

    Uses the Wallace rule for short oligonucleotides
    and the nearest-neighbor method for longer
    oligonucleotides.

    For sequences longer than 60 bases, the function
    returns 0.0 because a single whole-sequence Tm is
    not an appropriate interpretation for long genomic
    sequences.
    """

    if not sequence:
        return 0.0

    if not validate_dna(sequence):
        raise ValueError("Invalid DNA sequence.")

    sequence = sequence.upper()

    if len(sequence) <= 13:

        a = sequence.count("A")
        t = sequence.count("T")
        g = sequence.count("G")
        c = sequence.count("C")

        tm = 2 * (a + t) + 4 * (g + c)

    elif len(sequence) <= 60:

        tm = Tm_NN(sequence)

    else:

        return 0.0

    return round(float(tm), 2)