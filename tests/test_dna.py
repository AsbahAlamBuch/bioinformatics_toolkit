import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.dna import (
    get_sequence_length,
    calculate_gc_content,
    reverse_complement,
    validate_dna, kmer_count

)


def test_sequence_length():
    assert get_sequence_length("ATGC") == 4


def test_gc_content():
    assert calculate_gc_content("GGCC") == 100


def test_reverse_complement():
    assert reverse_complement("ATGC") == "GCAT"


def test_valid_dna():
    assert validate_dna("ATGC") is True


def test_invalid_dna():
    assert validate_dna("ATGX") is False

def test_kmer_count():
    expected = {
        "AT": 1,
        "TG": 2,
        "GT": 1,
    }

    assert kmer_count("ATGTG", 2) == expected
