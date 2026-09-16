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
    kmer_count,
    find_restriction_sites,
    calculate_melting_temperature
)

from src.protein import (
    translate_rna,
    amino_acid_composition,
    molecular_weight,
    codon_usage,
    calculate_gravy,
    calculate_pI,
    get_protein_length,
    calculate_instability_index,
    validate_protein
)

from src.fasta import read_fasta

from src.alignment import (
    global_alignment,
    local_alignment,
    hamming_distance,
    alignment_statistics
)


def analyze_sequence(
    dna: str,
    motif: str = "GTA",
    k: int = 3
) -> dict:
    """
    Analyze a DNA sequence and return biological properties.

    Parameters
    ----------
    dna : str
        DNA sequence to analyze.

    motif : str, optional
        DNA motif to search for.
        Default is "GTA".

    k : int, optional
        K-mer size.
        Default is 3.
    """

    dna = dna.upper()

    if not validate_dna(dna):
        raise ValueError("Invalid DNA sequence.")

    if not motif:
        raise ValueError("Motif cannot be empty.")

    motif = motif.upper()

    if not validate_dna(motif):
        raise ValueError("Invalid DNA motif.")

    if k <= 0:
        raise ValueError("k must be greater than 0.")

    rna = transcribe_dna(dna)

    translated_protein = translate_rna(rna)

    # Keep the complete translated sequence, including stop codons.
    protein = translated_protein.split("*")[0]

    if not validate_protein(protein):
        raise ValueError("Invalid protein sequence.")

    return {
        "sequence_length": get_sequence_length(dna),

        "gc_content": round(
            calculate_gc_content(dna),
            2
        ),

        "melting_temperature": calculate_melting_temperature(dna),

        "nucleotide_counts": count_nucleotides(dna),

        "reverse_sequence": reverse_sequence(dna),

        "reverse_complement": reverse_complement(dna),

        "rna_sequence": rna,

        "protein_sequence": translated_protein,

        "protein_length": get_protein_length(protein),

        "molecular_weight": molecular_weight(protein),

        "gravy_score": calculate_gravy(protein),

        "isoelectric_point": calculate_pI(protein),

        "instability_index": calculate_instability_index(protein),

        "amino_acid_composition": amino_acid_composition(protein),

        "codon_usage": codon_usage(rna),

        "motif_search": {
            "motif": motif,
            "positions": find_motif(dna, motif)
        },

        "kmer_frequencies": kmer_count(dna, k),

        "orfs": find_orfs(dna),

        "restriction_sites": find_restriction_sites(dna)
    }


def main():

    records = read_fasta("data/sample.fasta")

    for record in records:

        dna = str(record.seq).upper()

        print("=" * 60)
        print("Sequence ID:", record.id)
        print("DNA:", dna)

        try:
            results = analyze_sequence(dna)

        except ValueError as error:
            print(error)
            continue

        print("\nSequence Statistics")
        print("-" * 30)

        print(
            "Length:",
            results["sequence_length"]
        )

        print(
            "GC Content:",
            results["gc_content"],
            "%"
        )

        print(
            "Melting Temperature:",
            results["melting_temperature"],
            "°C"
        )

        print("\nNucleotide Counts")
        print("-" * 30)

        for nucleotide, count in (
            results["nucleotide_counts"].items()
        ):
            print(
                f"{nucleotide}: {count}"
            )

        print("\nReverse Sequence:")
        print(
            results["reverse_sequence"]
        )

        print("\nReverse Complement:")
        print(
            results["reverse_complement"]
        )

        print("\nRNA Sequence:")
        print(
            results["rna_sequence"]
        )

        print("\nProtein:")
        print(
            results["protein_sequence"]
        )

        print(
            "Protein Length:",
            results["protein_length"]
        )

        print("\nMolecular Weight:")
        print(
            results["molecular_weight"],
            "Da"
        )

        print("\nGRAVY Score:")
        print(
            results["gravy_score"]
        )

        print("\nEstimated pI:")
        print(
            results["isoelectric_point"]
        )

        print("\nInstability Index:")
        print(
            results["instability_index"]
        )

        print("\nAmino Acid Composition")
        print("-" * 30)

        for amino_acid, count in sorted(
            results["amino_acid_composition"].items()
        ):
            print(
                f"{amino_acid}: {count}"
            )

        print("\nCodon Usage")
        print("-" * 30)

        for codon, count in sorted(
            results["codon_usage"].items()
        ):
            print(
                f"{codon}: {count}"
            )

        print("\nMotif Search")
        print("-" * 30)

        motif_data = results["motif_search"]

        print(
            "Motif:",
            motif_data["motif"]
        )

        if motif_data["positions"]:

            print(
                "Found at positions:",
                motif_data["positions"]
            )

            print(
                "Number of matches:",
                len(motif_data["positions"])
            )

        else:

            print(
                "Motif not found"
            )

        print("\n3-mer Frequencies")
        print("-" * 30)

        for kmer, count in sorted(
            results["kmer_frequencies"].items()
        ):
            print(
                f"{kmer}: {count}"
            )

        print("\nORFs")
        print("-" * 30)

        if results["orfs"]:

            for start, end, sequence in (
                results["orfs"]
            ):

                print(
                    f"Start: {start}, End: {end}"
                )

                print(sequence)

        else:

            print(
                "No ORFs found"
            )

        print("\nRestriction Enzyme Sites")
        print("-" * 30)

        if results["restriction_sites"]:

            for enzyme, positions in (
                results["restriction_sites"].items()
            ):

                print(
                    f"{enzyme}: {positions}"
                )

        else:

            print(
                "No restriction sites found"
            )

    print("\n" + "=" * 60)
    print("GLOBAL ALIGNMENT")
    print("=" * 60)

    seq1 = "ATGCT"
    seq2 = "ATGTT"

    aligned_seq1, aligned_seq2, score = (
        global_alignment(seq1, seq2)
    )

    print(aligned_seq1)
    print(aligned_seq2)
    print("Score:", score)

    stats = alignment_statistics(
        aligned_seq1,
        aligned_seq2
    )

    print("\nAlignment Statistics")
    print("-" * 30)

    for key, value in stats.items():

        print(
            f"{key}: {value}"
        )

    print("\n" + "=" * 60)
    print("LOCAL ALIGNMENT")
    print("=" * 60)

    aligned_seq1, aligned_seq2, score = (
        local_alignment(seq1, seq2)
    )

    print(aligned_seq1)
    print(aligned_seq2)
    print("Score:", score)

    print("\n" + "=" * 60)
    print("HAMMING DISTANCE")
    print("=" * 60)

    hamming_seq1 = "ATGC"
    hamming_seq2 = "ATGT"

    distance = hamming_distance(
        hamming_seq1,
        hamming_seq2
    )

    print(
        "Sequence 1:",
        hamming_seq1
    )

    print(
        "Sequence 2:",
        hamming_seq2
    )

    print(
        "Distance:",
        distance
    )


if __name__ == "__main__":
    main()