import pytest

from src.dna import (
    get_sequence_length,
    validate_dna,
    calculate_gc_content,
    count_nucleotides,
    reverse_sequence,
    reverse_complement,
    transcribe_dna,
    find_motif,
    find_orfs,
    kmer_count,
    find_restriction_sites,
    calculate_melting_temperature,
)


def test_get_sequence_length():
    assert get_sequence_length("ATGC") == 4


def test_validate_dna():
    assert validate_dna("ATGC") is True
    assert validate_dna("ATGCXYZ") is False
    assert validate_dna("") is False


def test_validate_dna_lowercase():
    assert validate_dna("atgc") is True


def test_calculate_gc_content():
    assert calculate_gc_content("ATGC") == 50.0


def test_calculate_gc_content_invalid():
    with pytest.raises(ValueError):
        calculate_gc_content("ATGX")


def test_count_nucleotides():
    assert count_nucleotides("AATGCC") == {
        "A": 2,
        "T": 1,
        "G": 1,
        "C": 2,
    }


def test_reverse_sequence():
    assert reverse_sequence("ATGC") == "CGTA"


def test_reverse_complement():
    assert reverse_complement("ATGC") == "GCAT"


def test_transcribe_dna():
    assert transcribe_dna("ATGC") == "AUGC"


def test_find_motif():
    assert find_motif("ATGATG", "ATG") == [0, 3]


def test_find_motif_overlapping():
    assert find_motif("AAAA", "AA") == [0, 1, 2]


def test_find_motif_invalid():
    with pytest.raises(ValueError):
        find_motif("ATGC", "XYZ")


def test_find_orfs():
    sequence = "ATGAAATAA"

    assert find_orfs(sequence) == [
        (0, 9, "ATGAAATAA")
    ]


def test_find_orfs_forward_strand():
    sequence = "ATGAAATAA"

    assert (
        0,
        9,
        "ATGAAATAA"
    ) in find_orfs(sequence)


def test_find_orfs_reverse_complement():
    sequence = "TTATTTCAT"

    orfs = find_orfs(sequence)

    assert (
        9,
        0,
        "ATGAAATAA"
    ) in orfs


def test_kmer_count():
    assert kmer_count("ATAT", 2) == {
        "AT": 2,
        "TA": 1,
    }


def test_kmer_count_single():
    assert kmer_count("ATGC", 1) == {
        "A": 1,
        "T": 1,
        "G": 1,
        "C": 1,
    }


def test_kmer_count_larger_than_sequence():
    assert kmer_count("ATGC", 5) == {}


def test_kmer_count_invalid_k():
    with pytest.raises(ValueError):
        kmer_count("ATGC", 0)


def test_find_restriction_sites():
    sequence = "GAATTC"

    assert find_restriction_sites(sequence) == {
        "EcoRI": [0]
    }


def test_find_restriction_sites_none():
    assert find_restriction_sites("ATGCATGC") == {}


def test_calculate_melting_temperature_short():
    assert calculate_melting_temperature("ATGC") == 12.0


def test_calculate_melting_temperature_empty():
    assert calculate_melting_temperature("") == 0.0


def test_calculate_melting_temperature_invalid():
    with pytest.raises(ValueError):
        calculate_melting_temperature("ATGX")