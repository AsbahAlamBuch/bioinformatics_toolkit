import pytest

from src.fasta import read_fasta


def test_empty_fasta(tmp_path):
    fasta_file = tmp_path / "empty.fasta"
    fasta_file.write_text("")

    with pytest.raises(ValueError):
        read_fasta(fasta_file)