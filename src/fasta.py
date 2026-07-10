from Bio import SeqIO


def read_fasta(filename):
    records = list(SeqIO.parse(filename, "fasta"))
    return records