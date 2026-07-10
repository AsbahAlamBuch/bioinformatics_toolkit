from Bio.Seq import Seq
from src.dna import *
from src.protein import*


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

