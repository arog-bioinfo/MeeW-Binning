"""Shared contracts and helpers for MeeW-Binning job scripts."""

import csv
import re
from pathlib import Path

FASTA_SUFFIXES = {".fa", ".fna", ".fasta"}
ASSIGNMENT_COLUMNS = [
    "tool",
    "source_bin",
    "source_path",
    "normalized_bin_id",
    "assigned_domain",
    "prok_bp",
    "euk_bp",
    "labeled_bp",
    "unlabeled_bp",
    "total_bp",
    "n_contigs",
    "n_labeled_contigs",
    "fallback_domain",
]

USER_IDENTIFIER_PATTERN = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]{0,127}")
GENERATED_MAG_PATTERN = re.compile(
    r"(?P<sample>[A-Za-z0-9][A-Za-z0-9._-]{0,127})_"
    r"(?P<domain>prok|euk|unknown)_bin_(?P<number>[0-9]{3,})"
)
GENERATED_MAG_MAX_LENGTH = 192


def validate_user_identifier(value: str) -> str:
    """Return an unchanged valid user identifier or raise ValueError."""
    if not isinstance(value, str) or not value.isascii():
        raise ValueError(f"Invalid user identifier: {value!r}")
    if USER_IDENTIFIER_PATTERN.fullmatch(value) is None:
        raise ValueError(f"Invalid user identifier: {value!r}")
    return value


def validate_generated_mag_identifier(value: str) -> str:
    """Return an unchanged valid generated MAG identifier or raise ValueError."""
    if not isinstance(value, str) or not value.isascii():
        raise ValueError(f"Invalid generated MAG identifier: {value!r}")
    if len(value) > GENERATED_MAG_MAX_LENGTH:
        raise ValueError(
            "Generated MAG identifier exceeds 192 characters: " f"{value!r}"
        )
    match = GENERATED_MAG_PATTERN.fullmatch(value)
    if match is None:
        raise ValueError(f"Invalid generated MAG identifier: {value!r}")
    validate_user_identifier(match.group("sample"))
    return value


def generated_mag_identifier(sample: str, domain: str, number: int) -> str:
    """Form and validate ``<sample>_(prok|euk|unknown)_bin_<NNN>`` unchanged."""
    validate_user_identifier(sample)
    if domain not in {"prok", "euk", "unknown"}:
        raise ValueError(f"Invalid generated MAG domain: {domain!r}")
    if not isinstance(number, int) or isinstance(number, bool) or number < 1:
        raise ValueError(f"Invalid generated MAG number: {number!r}")
    return validate_generated_mag_identifier(f"{sample}_{domain}_bin_{number:03d}")


def normalize_domain(stage1: str, stage2: str, fallback_domain: str) -> str:
    """Map Tiara labels to prok/euk with a configured fallback."""
    text = f"{stage1} {stage2}".strip().lower()
    compact = text.replace("_", " ").replace("-", " ")

    prok_aliases = (
        "prokaryote",
        "prokarya",
        "bacteria",
        "archaea",
        "archae",
        "prok",
    )
    euk_aliases = (
        "eukaryote",
        "eukarya",
        "euk",
        "organelle",
        "mitochond",
        "plastid",
    )

    if any(alias in compact for alias in prok_aliases):
        return "prok"
    if any(alias in compact for alias in euk_aliases):
        return "euk"
    return fallback_domain


def read_tiara_labels(tsv_path: Path, fallback_domain: str) -> dict[str, str]:
    def normalize_header_cell(value: str) -> str:
        return " ".join(value.strip().lower().replace("_", " ").split())

    def is_tiara_header(row: list[str]) -> bool:
        normalized = [normalize_header_cell(cell) for cell in row]
        if not normalized:
            return False

        first = normalized[0]
        # Tiara output headers include wording like "sequence id" and stage columns.
        known_first = {"sequence id", "sequence", "contig", "contig id", "name"}
        if first in known_first:
            return True

        return any(
            cell in {"class", "class 1", "class 2", "stage 1", "stage 2"}
            for cell in normalized[1:3]
        )

    labels: dict[str, str] = {}
    with tsv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row:
                continue
            if is_tiara_header(row):
                continue

            contig_id = row[0].strip()
            stage1 = row[1].strip() if len(row) > 1 else ""
            stage2 = row[2].strip() if len(row) > 2 else ""
            if not contig_id:
                continue
            labels[contig_id] = normalize_domain(stage1, stage2, fallback_domain)
    return labels


def select_domain(bp_by_domain: dict[str, int], fallback_domain: str) -> str:
    prok_bp = bp_by_domain.get("prok", 0)
    euk_bp = bp_by_domain.get("euk", 0)
    if prok_bp > euk_bp:
        return "prok"
    if euk_bp > prok_bp:
        return "euk"
    return fallback_domain


def unique_fasta_lengths(fasta_path: Path) -> dict[str, int]:
    """Read FASTA lengths while rejecting duplicate or missing identifiers."""
    lengths: dict[str, int] = {}
    current_id: str | None = None
    current_len = 0

    with fasta_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    lengths[current_id] = current_len
                header = line[1:].split()
                current_id = header[0] if header else ""
                if not current_id:
                    raise ValueError(f"Missing FASTA identifier in {fasta_path}.")
                if current_id in lengths:
                    raise ValueError(
                        f"Duplicate FASTA identifier '{current_id}' in {fasta_path}."
                    )
                current_len = 0
            else:
                if current_id is None:
                    raise ValueError(
                        f"Sequence data appears before the first header in {fasta_path}."
                    )
                current_len += len(line)

    if current_id is not None:
        if current_id in lengths:
            raise ValueError(
                f"Duplicate FASTA identifier '{current_id}' in {fasta_path}."
            )
        lengths[current_id] = current_len
    return lengths


def serialize_provenance_path(path: Path, workdir: Path) -> str:
    """Serialize a canonical path relatively when it is inside the workdir."""
    resolved_path = path.resolve()
    resolved_workdir = workdir.resolve()
    try:
        return resolved_path.relative_to(resolved_workdir).as_posix()
    except ValueError:
        return resolved_path.as_posix()


def direct_child_fastas(directory: Path) -> list[Path]:
    """Return direct-child FASTAs in deterministic filename order."""
    return sorted(
        (
            path
            for path in directory.iterdir()
            if path.is_file() and path.suffix.lower() in FASTA_SUFFIXES
        ),
        key=lambda path: path.name,
    )
