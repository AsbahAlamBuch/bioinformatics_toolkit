from Bio import SeqIO


def read_fasta(filename):
    """
    Read sequences from a FASTA file.

    Returns a list of Biopython SeqRecord objects.
    """

    with open(filename, "r") as fasta_file:
        records = list(
            SeqIO.parse(
                fasta_file,
                "fasta"
            )
        )

    if not records:
        raise ValueError(
            "FASTA file contains no sequences."
        )

    return records