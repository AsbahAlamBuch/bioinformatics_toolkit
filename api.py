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
    description="A modular toolkit for DNA and protein sequence analysis.",
    version="1.0.0"
)


# Maximum FASTA upload size: 10 MB
MAX_FASTA_SIZE = 10 * 1024 * 1024


class SequenceRequest(BaseModel):
    sequence: str
    motif: str = Field(default="GTA", min_length=1)
    k: int = Field(default=3, gt=0)


@app.get("/")
def root():
    return FileResponse(
        Path(__file__).parent / "frontend" / "index.html"
    )


@app.post("/analyze")
def analyze(request: SequenceRequest):

    try:

        results = analyze_sequence(
            request.sequence,
            motif=request.motif,
            k=request.k
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

    if not motif:

        raise HTTPException(
            status_code=400,
            detail="Motif cannot be empty."
        )

    motif = motif.upper()

    if not all(
        base in "ATGC"
        for base in motif
    ):

        raise HTTPException(
            status_code=400,
            detail="Invalid DNA motif."
        )

    if k <= 0:

        raise HTTPException(
            status_code=400,
            detail="k must be greater than 0."
        )

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
                        detail="FASTA file is too large. Maximum size is 10 MB."
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