from Bio.Seq import Seq

def get_sequence_length(sequence):
    return len(sequence)

def count_nucleotides(sequence):
    print("A: ", sequence.count("A"))
    print("T: ", sequence.count("T"))
    print("G: ", sequence.count("G"))
    print("C: ", sequence.count("C"))

#main function start

def translate_rna(rna):
    return rna.translate()


def transcribe_dna(sequence):
    rna = sequence.replace("T", "U")
    return rna


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


def reverse_sequence(sequence):
    return sequence [::-1]
def calculate_gc_content(sequence):
    gc = sequence.count("G") + sequence.count("C")
    return(gc/len(sequence))*100

dna= Seq("ATGCGTAGCTA")

length= get_sequence_length(dna)
gc_content = calculate_gc_content(dna)

rna=transcribe_dna(dna)

protein= translate_rna(rna)

reversed_dna= reverse_sequence(dna)

reverse_comp = reverse_complement(dna)

print("DNA Sequence: ", dna)
print("Length: ", length)
print("Reversed Sequence: ", reversed_dna)
print("GC Content: ", round(gc_content, 2), "%")
print("Reverse Complement : ", reverse_comp)
print("RNA Sequence: ", rna)
print("Protein :", protein)


count_nucleotides(dna)

