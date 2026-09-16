from fastapi.testclient import TestClient

from api import app


client = TestClient(app)


def test_root():
    response = client.get("/")

    assert response.status_code == 200

    assert "Bioinformatics Toolkit" in response.text


def test_analyze_valid_sequence():
    response = client.post(
        "/analyze",
        json={
            "sequence": "ATGCGTAGCTAGCTAGCTAGCTAGCTAGCTA"
        }
    )

    assert response.status_code == 200

    data = response.json()

    assert "sequence_length" in data
    assert "gc_content" in data
    assert "rna_sequence" in data
    assert "protein_sequence" in data
    assert "molecular_weight" in data
    assert "gravy_score" in data
    assert "isoelectric_point" in data
    assert "instability_index" in data
    assert "orfs" in data
    assert "restriction_sites" in data


def test_analyze_custom_motif():
    response = client.post(
        "/analyze",
        json={
            "sequence": "ATGCGTAGCTAG",
            "motif": "ATG",
            "k": 3
        }
    )

    assert response.status_code == 200

    data = response.json()

    assert data["motif_search"]["motif"] == "ATG"
    assert data["motif_search"]["positions"] == [0]


def test_analyze_custom_kmer():
    response = client.post(
        "/analyze",
        json={
            "sequence": "ATGCGTAG",
            "motif": "GTA",
            "k": 4
        }
    )

    assert response.status_code == 200

    data = response.json()

    assert data["kmer_frequencies"] == {
        "ATGC": 1,
        "TGCG": 1,
        "GCGT": 1,
        "CGTA": 1,
        "GTAG": 1
    }


def test_analyze_invalid_motif():
    response = client.post(
        "/analyze",
        json={
            "sequence": "ATGCGTAG",
            "motif": "XYZ",
            "k": 3
        }
    )

    assert response.status_code == 400

    assert response.json() == {
        "detail": "Invalid DNA motif."
    }


def test_analyze_invalid_kmer():
    response = client.post(
        "/analyze",
        json={
            "sequence": "ATGCGTAG",
            "motif": "GTA",
            "k": 0
        }
    )

    assert response.status_code == 422


def test_analyze_invalid_sequence():
    response = client.post(
        "/analyze",
        json={
            "sequence": "ATGXYZ123"
        }
    )

    assert response.status_code == 400

    assert response.json() == {
        "detail": "Invalid DNA sequence."
    }


def test_analyze_empty_sequence():
    response = client.post(
        "/analyze",
        json={
            "sequence": ""
        }
    )

    assert response.status_code == 400


def test_analyze_lowercase_sequence():
    response = client.post(
        "/analyze",
        json={
            "sequence": "atgcgtagctag"
        }
    )

    assert response.status_code == 200

    data = response.json()

    assert data["sequence_length"] == 12
    assert data["rna_sequence"] == "AUGCGUAGCUAG"


def test_analyze_preserves_stop_codon():
    response = client.post(
        "/analyze",
        json={
            "sequence": "ATGAAATAA",
            "motif": "ATG",
            "k": 3
        }
    )

    assert response.status_code == 200

    data = response.json()

    assert data["protein_sequence"] == "MK*"
    assert data["protein_length"] == 2


def test_analyze_fasta():
    with open(
        "data/sample.fasta",
        "rb"
    ) as fasta_file:

        response = client.post(
            "/analyze-fasta",
            files={
                "file": (
                    "sample.fasta",
                    fasta_file,
                    "text/plain"
                )
            }
        )

    assert response.status_code == 200

    data = response.json()

    assert data["filename"] == "sample.fasta"
    assert data["motif"] == "GTA"
    assert data["k"] == 3
    assert data["sequence_count"] == 3

    assert len(data["sequences"]) == 3

    assert data["sequences"][0]["id"] == "Human_Gene"


def test_analyze_fasta_custom_parameters():
    with open(
        "data/sample.fasta",
        "rb"
    ) as fasta_file:

        response = client.post(
            "/analyze-fasta",
            files={
                "file": (
                    "sample.fasta",
                    fasta_file,
                    "text/plain"
                )
            },
            data={
                "motif": "ATG",
                "k": "4"
            }
        )

    assert response.status_code == 200

    data = response.json()

    assert data["motif"] == "ATG"
    assert data["k"] == 4
    assert data["sequence_count"] == 3

    for sequence in data["sequences"]:
        assert "analysis" in sequence

        analysis = sequence["analysis"]

        assert analysis["motif_search"]["motif"] == "ATG"

        assert all(
            len(kmer) == 4
            for kmer in analysis["kmer_frequencies"]
        )


def test_analyze_fasta_invalid_motif():
    with open(
        "data/sample.fasta",
        "rb"
    ) as fasta_file:

        response = client.post(
            "/analyze-fasta",
            files={
                "file": (
                    "sample.fasta",
                    fasta_file,
                    "text/plain"
                )
            },
            data={
                "motif": "XYZ",
                "k": "3"
            }
        )

    assert response.status_code == 400

    assert response.json() == {
        "detail": "Invalid DNA motif."
    }


def test_analyze_fasta_invalid_kmer():
    with open(
        "data/sample.fasta",
        "rb"
    ) as fasta_file:

        response = client.post(
            "/analyze-fasta",
            files={
                "file": (
                    "sample.fasta",
                    fasta_file,
                    "text/plain"
                )
            },
            data={
                "motif": "ATG",
                "k": "0"
            }
        )

    assert response.status_code == 400

    assert response.json() == {
        "detail": "k must be greater than 0."
    }


def test_analyze_fasta_invalid_extension():
    response = client.post(
        "/analyze-fasta",
        files={
            "file": (
                "sample.txt",
                b">test\nATGC",
                "text/plain"
            )
        }
    )

    assert response.status_code == 400

    assert response.json() == {
        "detail": "Please upload a FASTA file (.fasta, .fa, or .fna)."
    }


def test_analyze_fasta_empty_file():
    response = client.post(
        "/analyze-fasta",
        files={
            "file": (
                "empty.fasta",
                b"",
                "text/plain"
            )
        }
    )

    assert response.status_code == 400

    assert response.json() == {
        "detail": "Uploaded FASTA file is empty."
    }


def test_analyze_fasta_no_sequences():
    response = client.post(
        "/analyze-fasta",
        files={
            "file": (
                "empty.fasta",
                b"",
                "text/plain"
            )
        }
    )

    assert response.status_code == 400


def test_analyze_fasta_invalid_content():
    response = client.post(
        "/analyze-fasta",
        files={
            "file": (
                "invalid.fasta",
                b"this is not a FASTA file",
                "text/plain"
            )
        }
    )

    assert response.status_code == 400