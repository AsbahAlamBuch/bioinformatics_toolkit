from Bio.Seq import Seq

def get_sequence_length(sequence):
    return len(sequence)

def count_nucleotides(sequence):
    print("A: ", sequence.count("A"))
    print("T: ", sequence.count("T"))
    print("G: ", sequence.count("G"))
    print("C: ", sequence.count("C"))

dna= Seq("ATGCGTAGCTA")

length= get_sequence_length(dna)

print("DNA Sequence: ", dna)
print("Length: ", length)

count_nucleotides(dna)

