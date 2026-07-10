from Bio.Seq import Seq

def get_sequence_length(sequence):
    return len(sequence)


def count_nucleotides(sequence):
    print("A: ", sequence.count("A"))
    print("T: ", sequence.count("T"))
    print("G: ", sequence.count("G"))
    print("C: ", sequence.count("C"))

def calculate_gc_content(sequence):
    gc = sequence.count("G") + sequence.count("C")
    return(gc/len(sequence))*100


def reverse_sequence(sequence):
    return sequence [::-1]

def reverse_complement(sequence):

    complements = {"A":"T",
                   "T": "A",
                   "G": "C",
                   "C": "G"}

    reverse = sequence[::-1]
    result = ""

    for base in reverse:
        result += complements[base]

    return result

def transcribe_dna(sequence):
    rna = sequence.replace("T", "U")
    return rna

