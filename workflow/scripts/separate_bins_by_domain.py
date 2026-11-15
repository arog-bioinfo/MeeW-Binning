"""Assign source bins to prok/euk/unknown domains using Tiara labels."""

import csv
import shutil
import sys
from collections import defaultdict
from pathlib import Path

_SCRIPTS_DIR = (
    Path(snakemake.scriptdir).resolve()
    if "snakemake" in globals() and hasattr(snakemake, "scriptdir")
    else Path(__file__).resolve().parent
)
sys.path.insert(0, str(_SCRIPTS_DIR))

from meew_binning_common import (  # noqa: E402
    generated_mag_identifier,
    read_tiara_labels,
    select_domain,
)


def fasta_lengths(fasta_path: Path) -> dict[str, int]:
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
                current_id = line[1:].split()[0]
                current_len = 0
            else:
                current_len += len(line)

    if current_id is not None:
        lengths[current_id] = current_len
    return lengths


def relative_path(path: Path, start: Path) -> str:
    """Return a POSIX-style path relative to start."""
    return path.resolve().relative_to(start.resolve()).as_posix()


def main() -> None:
    bins_dir_value = getattr(snakemake.params, "bins_dir", None)
    if bins_dir_value:
        bins_dir = Path(bins_dir_value)
    else:
        bins_dir = Path(snakemake.input.bins)
    tiara_tsv = Path(snakemake.input.tiara_tsv)
    prok_dir = Path(snakemake.output.prok_bins)
    euk_dir = Path(snakemake.output.euk_bins)
    unknown_dir = Path(snakemake.output.unknown_bins)
    summary_path = Path(snakemake.output.summary)
    annotation_samples_path = Path(snakemake.output.annotation_samples)
    rename_map_path = Path(snakemake.output.rename_map)
    fallback_domain = str(snakemake.params.fallback_domain)
    sample_prefix = str(snakemake.wildcards.sample)
    workdir = Path.cwd()

    domain_dirs = {"prok": prok_dir, "euk": euk_dir, "unknown": unknown_dir}
    for domain_dir in domain_dirs.values():
        if domain_dir.exists():
            shutil.rmtree(domain_dir)
        domain_dir.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    annotation_samples_path.parent.mkdir(parents=True, exist_ok=True)
    rename_map_path.parent.mkdir(parents=True, exist_ok=True)

    tiara_labels = read_tiara_labels(tiara_tsv, fallback_domain)
    fasta_suffixes = {".fa", ".fna", ".fasta"}
    bin_fastas = sorted(
        path
        for path in bins_dir.iterdir()
        if path.is_file() and path.suffix.lower() in fasta_suffixes
    )

    rows: list[dict[str, str | int]] = []
    annotation_rows: list[dict[str, str]] = []
    rename_rows: list[dict[str, str]] = []
    domain_counts = defaultdict(int)
    for fasta_path in bin_fastas:
        lengths = fasta_lengths(fasta_path)
        bp_by_domain = defaultdict(int)
        labeled_bp = 0
        unlabeled_bp = 0

        for contig_id, contig_bp in lengths.items():
            domain = tiara_labels.get(contig_id)
            if domain is None:
                unlabeled_bp += contig_bp
                continue
            bp_by_domain[domain] += contig_bp
            labeled_bp += contig_bp

        assigned_domain = select_domain(bp_by_domain, fallback_domain)
        target_dir = domain_dirs[assigned_domain]
        domain_counts[assigned_domain] += 1
        sample_name = generated_mag_identifier(
            sample_prefix, assigned_domain, domain_counts[assigned_domain]
        )
        standardized_name = f"{sample_name}.fa"
        standardized_path = target_dir / standardized_name
        shutil.copy2(fasta_path, standardized_path)

        original_path = relative_path(fasta_path, workdir)
        standardized_relative_path = relative_path(standardized_path, workdir)

        total_bp = labeled_bp + unlabeled_bp
        rows.append(
            {
                "bin": fasta_path.name,
                "original_bin": fasta_path.name,
                "original_path": original_path,
                "standardized_bin": standardized_name,
                "standardized_path": standardized_relative_path,
                "assigned_domain": assigned_domain,
                "prok_bp": bp_by_domain.get("prok", 0),
                "euk_bp": bp_by_domain.get("euk", 0),
                "labeled_bp": labeled_bp,
                "unlabeled_bp": unlabeled_bp,
                "total_bp": total_bp,
                "n_contigs": len(lengths),
                "n_labeled_contigs": sum(1 for cid in lengths if cid in tiara_labels),
                "fallback_domain": fallback_domain,
            }
        )
        if assigned_domain != "unknown":
            annotation_rows.append(
                {
                    "sample": sample_name,
                    "path": standardized_relative_path,
                    "domain": assigned_domain,
                    "genome_type": "mag",
                }
            )
        rename_rows.append(
            {
                "original_bin": fasta_path.name,
                "original_path": original_path,
                "standardized_bin": standardized_name,
                "standardized_path": standardized_relative_path,
                "domain": assigned_domain,
            }
        )

    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = [
            "bin",
            "original_bin",
            "original_path",
            "standardized_bin",
            "standardized_path",
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
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    with annotation_samples_path.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = ["sample", "path", "domain", "genome_type"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in annotation_rows:
            writer.writerow(row)

    with rename_map_path.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = [
            "original_bin",
            "original_path",
            "standardized_bin",
            "standardized_path",
            "domain",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rename_rows:
            writer.writerow(row)


if __name__ == "__main__":
    main()
