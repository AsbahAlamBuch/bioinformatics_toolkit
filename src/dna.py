def get_sequence_length(sequence):
    return len(sequence)


def count_nucleotides(sequence):
    print("A:", sequence.count("A"))
    print("T:", sequence.count("T"))
    print("G:", sequence.count("G"))
    print("C:", sequence.count("C"))


def calculate_gc_content(sequence):
    gc = sequence.count("G") + sequence.count("C")
    return (gc / len(sequence)) * 100


def reverse_sequence(sequence):
    return sequence[::-1]


def reverse_complement(sequence):
    complements = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G"
    }

    reverse = sequence[::-1]
    result = ""

    for base in reverse:
        result += complements[base]

    return result


def transcribe_dna(sequence):
    return sequence.replace("T", "U")


def find_motif(sequence, motif):
    positions = []

    for i in range(len(sequence) - len(motif) + 1):
        if sequence[i:i + len(motif)] == motif:
            positions.append(i + 1)

    return positions

def validate_dna(sequence):
    valid_bases = {"A", "T", "G", "C"}

    for i, base in enumerate(sequence):
        if base not in valid_bases:
            print(f"Invalid nucleotide '{base}' at position {i + 1}")
            return False

    print("DNA sequence is valid.")
    return True

def find_orfs(sequence):
    pass

def find_orfs(sequence):
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]

    orfs = []

    for i in range(len(sequence) - 2):
        codon = sequence[i:i+3]

        if codon == start_codon:
            for j in range(i + 3, len(sequence) - 2, 3):
                stop = sequence[j:j+3]

                if stop in stop_codons:
                    orfs.append((i + 1, j + 3, sequence[i:j+3]))
                    break

    return orfs