from typing import List, Dict


def _normalize_sequences(seq1: str, seq2: str) -> tuple[str, str]:
    """Validate and normalize two sequences for alignment."""

    if not isinstance(seq1, str) or not isinstance(seq2, str):
        raise ValueError("Sequences must be strings.")

    return seq1.upper(), seq2.upper()


def global_alignment(
    seq1: str,
    seq2: str,
    match: int = 1,
    mismatch: int = -1,
    gap: int = -2,
):
    """
    Needleman-Wunsch Global Sequence Alignment.

    Empty sequences are allowed.
    """

    seq1, seq2 = _normalize_sequences(seq1, seq2)

    rows = len(seq1) + 1
    cols = len(seq2) + 1

    score_matrix: List[List[int]] = [
        [0 for _ in range(cols)]
        for _ in range(rows)
    ]

    for i in range(rows):
        score_matrix[i][0] = i * gap

    for j in range(cols):
        score_matrix[0][j] = j * gap

    for i in range(1, rows):
        for j in range(1, cols):

            if seq1[i - 1] == seq2[j - 1]:
                diagonal = score_matrix[i - 1][j - 1] + match
            else:
                diagonal = score_matrix[i - 1][j - 1] + mismatch

            up = score_matrix[i - 1][j] + gap
            left = score_matrix[i][j - 1] + gap

            score_matrix[i][j] = max(
                diagonal,
                up,
                left
            )

    aligned_seq1 = ""
    aligned_seq2 = ""

    i = rows - 1
    j = cols - 1

    while i > 0 or j > 0:

        if i > 0 and j > 0:

            if seq1[i - 1] == seq2[j - 1]:
                current_score = match
            else:
                current_score = mismatch

            if score_matrix[i][j] == (
                score_matrix[i - 1][j - 1] + current_score
            ):
                aligned_seq1 = seq1[i - 1] + aligned_seq1
                aligned_seq2 = seq2[j - 1] + aligned_seq2

                i -= 1
                j -= 1
                continue

        if i > 0 and score_matrix[i][j] == (
            score_matrix[i - 1][j] + gap
        ):
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2

            i -= 1

        else:
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2

            j -= 1

    alignment_score = score_matrix[rows - 1][cols - 1]

    return aligned_seq1, aligned_seq2, alignment_score


def local_alignment(
    seq1: str,
    seq2: str,
    match: int = 1,
    mismatch: int = -1,
    gap: int = -2,
):
    """
    Smith-Waterman Local Sequence Alignment.

    Empty sequences are allowed.

    Returns
    -------
    tuple
        (aligned_seq1, aligned_seq2, score)
    """

    seq1, seq2 = _normalize_sequences(seq1, seq2)

    rows = len(seq1) + 1
    cols = len(seq2) + 1

    score_matrix: List[List[int]] = [
        [0 for _ in range(cols)]
        for _ in range(rows)
    ]

    max_score = 0
    max_position = (0, 0)

    for i in range(1, rows):
        for j in range(1, cols):

            if seq1[i - 1] == seq2[j - 1]:
                diagonal = score_matrix[i - 1][j - 1] + match
            else:
                diagonal = score_matrix[i - 1][j - 1] + mismatch

            up = score_matrix[i - 1][j] + gap
            left = score_matrix[i][j - 1] + gap

            score_matrix[i][j] = max(
                0,
                diagonal,
                up,
                left
            )

            if score_matrix[i][j] > max_score:
                max_score = score_matrix[i][j]
                max_position = (i, j)

    aligned_seq1 = ""
    aligned_seq2 = ""

    i, j = max_position

    while i > 0 and j > 0 and score_matrix[i][j] > 0:

        if seq1[i - 1] == seq2[j - 1]:
            current_score = match
        else:
            current_score = mismatch

        if score_matrix[i][j] == (
            score_matrix[i - 1][j - 1] + current_score
        ):
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2

            i -= 1
            j -= 1

        elif score_matrix[i][j] == (
            score_matrix[i - 1][j] + gap
        ):
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2

            i -= 1

        else:
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2

            j -= 1

    return aligned_seq1, aligned_seq2, max_score


def hamming_distance(seq1: str, seq2: str) -> int:
    """
    Calculate the Hamming distance between two sequences.

    The sequences must be non-empty and have equal length.
    """

    if not isinstance(seq1, str) or not isinstance(seq2, str):
        raise ValueError("Sequences must be strings.")

    if not seq1 or not seq2:
        raise ValueError("Sequences cannot be empty.")

    seq1 = seq1.upper()
    seq2 = seq2.upper()

    if len(seq1) != len(seq2):
        raise ValueError(
            "Hamming distance requires sequences of equal length."
        )

    return sum(
        nucleotide1 != nucleotide2
        for nucleotide1, nucleotide2 in zip(seq1, seq2)
    )


def alignment_statistics(
    aligned_seq1: str,
    aligned_seq2: str,
) -> Dict[str, float]:
    """
    Calculate statistics for two aligned sequences.
    """

    if not isinstance(aligned_seq1, str) or not isinstance(
        aligned_seq2, str
    ):
        raise ValueError("Aligned sequences must be strings.")

    if len(aligned_seq1) != len(aligned_seq2):
        raise ValueError(
            "Aligned sequences must have equal length."
        )

    matches = 0
    mismatches = 0
    gaps = 0

    for base1, base2 in zip(
        aligned_seq1.upper(),
        aligned_seq2.upper()
    ):

        if base1 == "-" or base2 == "-":
            gaps += 1

        elif base1 == base2:
            matches += 1

        else:
            mismatches += 1

    alignment_length = len(aligned_seq1)

    if alignment_length > 0:
        identity_percent = (
            matches / alignment_length
        ) * 100
    else:
        identity_percent = 0.0

    return {
        "identity_percent": round(identity_percent, 2),
        "matches": matches,
        "mismatches": mismatches,
        "gaps": gaps,
        "alignment_length": alignment_length,
    }