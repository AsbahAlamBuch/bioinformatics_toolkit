import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.protein import (
    codon_usage,
    amino_acid_composition,
    molecular_weight,
    validate_protein,
    calculate_gravy,
    calculate_pI,
    get_protein_length,
    calculate_instability_index,
)

from src.protein import molecular_weight

def test_molecular_weight():
    assert molecular_weight("AG") == 146.14


def test_codon_usage():
    rna = "AUGGCUGCUAAA"

    expected = {
        "AUG": 1,
        "GCU": 2,
        "AAA": 1
    }

    assert codon_usage(rna) == expected

def test_amino_acid_composition():
    protein = "MAGAA"

    expected = {
        "M": 1,
        "A": 3,
        "G": 1
    }

    assert amino_acid_composition(protein) == expected

def test_validate_protein():
    assert validate_protein("MAGK") is True
    assert validate_protein("MAGX") is False


def test_gravy():
    assert calculate_gravy("AAAA") == 1.8


def test_protein_length():
    assert get_protein_length("MAGK") == 4


def test_pI():
    assert calculate_pI("K") == 8.75


def test_instability_index():
    assert calculate_instability_index("MRAIYISRARA") == 19.75

def test_invalid_instability_index():
    import pytest
    from src.protein import calculate_instability_index

    with pytest.raises(ValueError):
        calculate_instability_index("MAGX")

def test_invalid_pI():
    import pytest
    from src.protein import calculate_pI

    with pytest.raises(ValueError):
        calculate_pI("MAGX")

def test_invalid_rna_translation():
    import pytest
    from src.protein import translate_rna

    with pytest.raises(ValueError):
        translate_rna("AUGXCC")

def test_empty_rna_translation():
    from src.protein import translate_rna

    assert translate_rna("") == ""

def test_incomplete_rna_translation():
    from src.protein import translate_rna

    assert translate_rna("AUGGCUA") == "MA"

def test_empty_codon_usage():
    from src.protein import codon_usage

    assert codon_usage("") == {}

def test_incomplete_codon_usage():
    from src.protein import codon_usage

    assert codon_usage("AUGGCUA") == {
        "AUG": 1,
        "GCU": 1,
    }

def test_amino_acid_composition_lowercase():
    assert amino_acid_composition("magaa") == {
        "M": 1,
        "A": 3,
        "G": 1
    }


def test_invalid_amino_acid_composition():
    import pytest

    with pytest.raises(ValueError):
        amino_acid_composition("MAGX")

def test_validate_empty_protein():
    assert validate_protein("") is False