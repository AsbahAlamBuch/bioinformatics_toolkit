from src.dna import (
    get_sequence_length,
    calculate_gc_content,
    reverse_sequence,
    reverse_complement,
    transcribe_dna,
    count_nucleotides,
    find_motif,
    validate_dna,
    find_orfs,
)

from src.protein import translate_rna
from src.fasta import read_fasta

from src.alignment import global_alignment

records = read_fasta("data/sample.fasta")

for record in records:
    dna = str(record.seq)
    print(record.id)
    print(dna)

    if not validate_dna(dna):
        continue



    length = get_sequence_length(dna)
    gc_content = calculate_gc_content(dna)

    rna = transcribe_dna(dna)
    protein = translate_rna(rna)

    reversed_dna = reverse_sequence(dna)
    reverse_comp = reverse_complement(dna)

    motif = "GTA"
    positions = find_motif(dna, motif)

    orfs = find_orfs(dna)

    print("=" * 60)
    print("Sequence ID:", record.id)
    print("DNA:", dna)
    print("Length:", length)
    print("GC Content:", round(gc_content, 2), "%")
    print("Reverse Sequence:", reversed_dna)
    print("Reverse Complement:", reverse_comp)
    print("RNA Sequence:", rna)
    print("Protein:", protein)

    print("Motif:", motif)

    if positions:
        print("Found at positions:", positions)
        print("Number of matches:", len(positions))
    else:
        print("Motif not found")

    count_nucleotides(dna)
    print()

    print("ORFs:")

    if orfs:
        for start, end, seq in orfs:
            print(f"Start: {start}, End: {end}")
            print(seq)
    else:
        print("No ORFs found")


seq1 = "ATGCT"
seq2 = "ATGTT"

aligned1, aligned2, score = global_alignment(seq1, seq2)

print(aligned1)
print(aligned2)
print(score)