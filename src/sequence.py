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

dna= Seq("ATGCGTAGCTA")

length= get_sequence_length(dna)
gc_content = calculate_gc_content(dna)

print("DNA Sequence: ", dna)
print("Length: ", length)

print("GC Content: ", round(gc_content, 2), "%")


count_nucleotides(dna)

