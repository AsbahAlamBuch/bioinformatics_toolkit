import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.protein import (
    codon_usage,
    amino_acid_composition,
)

from src.protein import molecular_weight

def test_molecular_weight():
    assert molecular_weight("AG") == 164.16


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
