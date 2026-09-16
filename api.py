from fastapi import FastAPI, HTTPException, UploadFile, File, Form
from fastapi.responses import FileResponse
from pydantic import BaseModel, Field
from pathlib import Path
import os
import tempfile

from main import analyze_sequence
from src.fasta import read_fasta


app = FastAPI(
    title="Bioinformatics Toolkit",
    description="A complete toolkit for DNA and protein sequence analysis.",
    version="1.0.0"
)


# -------------------------------------------------------------------
# API limits
# -------------------------------------------------------------------

MAX_DNA_LENGTH = 100_000
MAX_MOTIF_LENGTH = 100
MAX_KMER_SIZE = 20
MAX_FASTA_SIZE = 10 * 1024 * 1024


# -------------------------------------------------------------------
# Request models
# -------------------------------------------------------------------

class SequenceRequest(BaseModel):
    sequence: str = Field(
        ...,
        max_length=MAX_DNA_LENGTH,
        description="DNA sequence containing only A, T, G and C."
    )

    motif: str = Field(
        default="GTA",
        max_length=MAX_MOTIF_LENGTH,
        description="DNA motif to search for."
    )

    k: int = Field(
        default=3,
        gt=0,
        le=MAX_KMER_SIZE,
        description="K-mer size."
    )


# -------------------------------------------------------------------
# Helper validation
# -------------------------------------------------------------------

def validate_motif(motif: str) -> str:
    """Validate and normalize a DNA motif."""

    if not motif:
        raise HTTPException(
            status_code=400,
            detail="Motif cannot be empty."
        )

    motif = motif.upper()

    if len(motif) > MAX_MOTIF_LENGTH:
        raise HTTPException(
            status_code=400,
            detail="Invalid DNA motif."
        )

    if not all(base in "ATGC" for base in motif):
        raise HTTPException(
            status_code=400,
            detail="Invalid DNA motif."
        )

    return motif


def validate_k(k: int) -> int:
    """Validate the k-mer size."""

    if k <= 0:
        raise HTTPException(
            status_code=400,
            detail="k must be greater than 0."
        )

    if k > MAX_KMER_SIZE:
        raise HTTPException(
            status_code=400,
            detail=f"k is too large. Maximum k is {MAX_KMER_SIZE}."
        )

    return k


# -------------------------------------------------------------------
# Routes
# -------------------------------------------------------------------

@app.get("/")
def root():
    """Serve the web frontend."""

    return FileResponse(
        Path(__file__).parent / "frontend" / "index.html"
    )


@app.get("/health")
def health():
    """Return a simple application health check."""

    return {
        "status": "healthy",
        "service": "Bioinformatics Toolkit",
        "version": "1.0.0"
    }


@app.post("/analyze")
def analyze(request: SequenceRequest):
    """Analyze a single DNA sequence."""

    if not request.sequence:
        raise HTTPException(
            status_code=400,
            detail="Invalid DNA sequence."
        )

    motif = validate_motif(request.motif)
    k = validate_k(request.k)

    sequence = request.sequence.upper()

    if len(sequence) > MAX_DNA_LENGTH:
        raise HTTPException(
            status_code=413,
            detail=(
                f"DNA sequence is too long. "
                f"Maximum length is {MAX_DNA_LENGTH} bases."
            )
        )

    try:

        results = analyze_sequence(
            sequence,
            motif=motif,
            k=k
        )

        return results

    except ValueError as error:

        raise HTTPException(
            status_code=400,
            detail=str(error)
        )


@app.post("/analyze-fasta")
async def analyze_fasta(
    file: UploadFile = File(...),
    motif: str = Form(default="GTA"),
    k: int = Form(default=3)
):
    """Analyze all DNA sequences contained in a FASTA file."""

    if not file.filename:

        raise HTTPException(
            status_code=400,
            detail="A FASTA filename is required."
        )

    if not file.filename.lower().endswith(
        (".fasta", ".fa", ".fna")
    ):

        raise HTTPException(
            status_code=400,
            detail="Please upload a FASTA file (.fasta, .fa, or .fna)."
        )

    motif = validate_motif(motif)
    k = validate_k(k)

    temporary_path = None

    try:

        file_descriptor, temporary_filename = tempfile.mkstemp(
            suffix=".fasta"
        )

        temporary_path = Path(
            temporary_filename
        )

        os.close(file_descriptor)

        total_size = 0

        with open(
            temporary_path,
            "wb"
        ) as temporary_file:

            while True:

                chunk = await file.read(1024 * 1024)

                if not chunk:
                    break

                total_size += len(chunk)

                if total_size > MAX_FASTA_SIZE:

                    raise HTTPException(
                        status_code=413,
                        detail=(
                            "FASTA file is too large. "
                            "Maximum size is 10 MB."
                        )
                    )

                temporary_file.write(chunk)

        if total_size == 0:

            raise HTTPException(
                status_code=400,
                detail="Uploaded FASTA file is empty."
            )

        records = read_fasta(
            temporary_path
        )

        results = []

        for record in records:

            dna = str(
                record.seq
            ).upper()

            if len(dna) > MAX_DNA_LENGTH:

                results.append({
                    "id": record.id,
                    "description": record.description,
                    "error": (
                        f"Sequence is too long. "
                        f"Maximum length is {MAX_DNA_LENGTH} bases."
                    )
                })

                continue

            try:

                analysis = analyze_sequence(
                    dna,
                    motif=motif,
                    k=k
                )

                results.append({
                    "id": record.id,
                    "description": record.description,
                    "analysis": analysis
                })

            except ValueError as error:

                results.append({
                    "id": record.id,
                    "description": record.description,
                    "error": str(error)
                })

        return {
            "filename": file.filename,
            "motif": motif,
            "k": k,
            "sequence_count": len(results),
            "sequences": results
        }

    except HTTPException:
        raise

    except ValueError as error:

        raise HTTPException(
            status_code=400,
            detail=str(error)
        )

    finally:

        if (
            temporary_path is not None
            and temporary_path.exists()
        ):

            try:

                os.remove(
                    temporary_path
                )

            except PermissionError:

                pass