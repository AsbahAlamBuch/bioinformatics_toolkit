import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.alignment import (
    global_alignment,
    local_alignment,
    hamming_distance,
    alignment_statistics,
)


def test_global_alignment():
    seq1, seq2, score = global_alignment("ATGCT", "ATGTT")

    assert seq1 == "ATGCT"
    assert seq2 == "ATGTT"
    assert score == 3


def test_local_alignment():
    seq1, seq2, score = local_alignment("ATGCT", "ATGTT")

    assert seq1 == "ATG"
    assert seq2 == "ATG"
    assert score == 3


def test_hamming_distance():
    assert hamming_distance("ATGC", "ATGT") == 1


def test_alignment_statistics():
    stats = alignment_statistics("ATGCT", "ATGTT")

    assert stats["identity_percent"] == 80.0
    assert stats["matches"] == 4
    assert stats["mismatches"] == 1
    assert stats["gaps"] == 0
    assert stats["alignment_length"] == 5


def test_hamming_distance_requires_equal_length():
    import pytest

    with pytest.raises(ValueError):
        hamming_distance("ATGC", "ATG")


def test_global_alignment_empty_sequences():
    seq1, seq2, score = global_alignment("", "")

    assert seq1 == ""
    assert seq2 == ""
    assert score == 0


def test_local_alignment_empty_sequences():
    seq1, seq2, score = local_alignment("", "")

    assert seq1 == ""
    assert seq2 == ""
    assert score == 0


def test_alignment_statistics_requires_equal_length():
    import pytest

    with pytest.raises(ValueError):
        alignment_statistics("ATGC", "ATG")


def test_alignment_statistics_empty():
    stats = alignment_statistics("", "")

    assert stats["identity_percent"] == 0.0
    assert stats["matches"] == 0
    assert stats["mismatches"] == 0
    assert stats["gaps"] == 0
    assert stats["alignment_length"] == 0


def test_global_alignment_empty_first_sequence():
    seq1, seq2, score = global_alignment("", "ATGC")

    assert seq1 == "----"
    assert seq2 == "ATGC"
    assert score == -8


def test_global_alignment_empty_second_sequence():
    seq1, seq2, score = global_alignment("ATGC", "")

    assert seq1 == "ATGC"
    assert seq2 == "----"
    assert score == -8


def test_local_alignment_no_matching_region():
    seq1, seq2, score = local_alignment("AAAA", "TTTT")

    assert seq1 == ""
    assert seq2 == ""
    assert score == 0


def test_hamming_distance_identical_sequences():
    assert hamming_distance("ATGC", "ATGC") == 0


def test_hamming_distance_case_insensitive():
    assert hamming_distance("atgc", "ATGT") == 1


def test_hamming_distance_empty_sequence():
    import pytest

    with pytest.raises(ValueError):
        hamming_distance("", "ATGC")


def test_alignment_statistics_case_insensitive():
    stats = alignment_statistics(
        "atgc",
        "ATGT"
    )

    assert stats["matches"] == 3
    assert stats["mismatches"] == 1
    assert stats["gaps"] == 0
    assert stats["alignment_length"] == 4
    assert stats["identity_percent"] == 75.0


def test_alignment_statistics_non_string():
    import pytest

    with pytest.raises(ValueError):
        alignment_statistics(
            "ATGC",
            None
        )