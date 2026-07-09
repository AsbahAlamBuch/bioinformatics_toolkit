from Bio.Seq import Seq

dna= Seq ("ATGCGTAGCTA")

print("DNA Sequence : ", dna)
print("Sequence Length: ", len(dna))
print("A :", dna.count("A"))
print("T :", dna.count("T"))
print("G :", dna.count("G"))
print("C :", dna.count("C"))