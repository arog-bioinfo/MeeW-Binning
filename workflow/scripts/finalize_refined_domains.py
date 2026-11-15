"""Materialize final prokaryotic, eukaryotic, and unknown bin directories."""

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
    ASSIGNMENT_COLUMNS,
    direct_child_fastas,
    generated_mag_identifier,
    read_tiara_labels,
    select_domain,
    serialize_provenance_path,
    unique_fasta_lengths,
)

INTEGER_COLUMNS = [
    "prok_bp",
    "euk_bp",
    "labeled_bp",
    "unlabeled_bp",
    "total_bp",
    "n_contigs",
    "n_labeled_contigs",
]


def parse_provenance_path(value: str, workdir: Path) -> Path:
    """Parse a relative-or-absolute provenance path to its canonical location."""
    if not value:
        raise ValueError("Provenance paths must be non-empty.")
    path = Path(value)
    if not path.is_absolute():
        path = workdir / path
    return path.resolve()


def summarize_bin(
    lengths: dict[str, int], tiara_labels: dict[str, str]
) -> dict[str, int]:
    """Calculate the shared Tiara base-pair and contig statistics."""
    bp_by_domain = defaultdict(int)
    labeled_bp = 0
    unlabeled_bp = 0
    for contig_id, contig_bp in lengths.items():
        tiara_domain = tiara_labels.get(contig_id)
        if tiara_domain is None:
            unlabeled_bp += contig_bp
            continue
        bp_by_domain[tiara_domain] += contig_bp
        labeled_bp += contig_bp
    return {
        "prok_bp": bp_by_domain.get("prok", 0),
        "euk_bp": bp_by_domain.get("euk", 0),
        "labeled_bp": labeled_bp,
        "unlabeled_bp": unlabeled_bp,
        "total_bp": labeled_bp + unlabeled_bp,
        "n_contigs": len(lengths),
        "n_labeled_contigs": sum(
            1 for contig_id in lengths if contig_id in tiara_labels
        ),
    }


def parse_nonnegative_integer(value: str, column: str, line_number: int) -> int:
    """Parse canonical non-negative integer assignment metadata."""
    if not value or not value.isascii() or not value.isdecimal():
        raise ValueError(
            f"Invalid non-negative integer '{value}' for {column} in bin-domain "
            f"assignment row {line_number}."
        )
    return int(value)


def validated_assignments(
    assignments_path: Path,
    tools: list[str],
    source_dirs: list[Path],
    sample: str,
    fallback_domain: str,
    tiara_labels: dict[str, str],
    workdir: Path,
) -> list[dict[str, str | int | Path]]:
    """Validate consolidated metadata against every declared source directory."""
    if len(tools) != len(source_dirs) or len(set(tools)) != len(tools):
        raise ValueError(
            "Selected binner tools and source directories must be unique, ordered, "
            "and have matching lengths."
        )

    expected_rows = []
    for tool, source_dir in zip(tools, source_dirs, strict=True):
        for index, fasta_path in enumerate(direct_child_fastas(source_dir), start=1):
            expected_rows.append((tool, source_dir, fasta_path, index))

    with assignments_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ASSIGNMENT_COLUMNS:
            raise ValueError(
                f"Invalid bin-domain assignment columns in {assignments_path}: "
                f"expected {ASSIGNMENT_COLUMNS}, observed {reader.fieldnames}."
            )
        raw_rows = list(reader)

    if len(raw_rows) != len(expected_rows):
        raise ValueError(
            f"Bin-domain assignment row count {len(raw_rows)} does not match the "
            f"{len(expected_rows)} direct-child source FASTAs."
        )

    validated_rows: list[dict[str, str | int | Path]] = []
    for line_number, (row, expected) in enumerate(
        zip(raw_rows, expected_rows, strict=True), start=2
    ):
        if None in row or any(value is None for value in row.values()):
            raise ValueError(
                f"Malformed bin-domain assignment row {line_number} in "
                f"{assignments_path}."
            )
        tool, source_dir, fasta_path, index = expected
        expected_source_path = serialize_provenance_path(fasta_path, workdir)
        expected_normalized_id = f"{sample}__{tool}__bin_{index:06d}"
        expected_identity = {
            "tool": tool,
            "source_bin": fasta_path.name,
            "source_path": expected_source_path,
            "normalized_bin_id": expected_normalized_id,
            "fallback_domain": fallback_domain,
        }
        for column, expected_value in expected_identity.items():
            if row[column] != expected_value:
                raise ValueError(
                    f"Inconsistent {column} in bin-domain assignment row "
                    f"{line_number}: expected {expected_value!r}, observed "
                    f"{row[column]!r}."
                )

        declared_source = parse_provenance_path(row["source_path"], workdir)
        if (
            declared_source.resolve().parent != source_dir.resolve()
            or declared_source.resolve() != fasta_path.resolve()
        ):
            raise ValueError(
                f"Bin-domain assignment source_path in row {line_number} is not "
                f"the expected direct child of binner '{tool}' directory "
                f"{source_dir}."
            )

        numeric_values = {
            column: parse_nonnegative_integer(row[column], column, line_number)
            for column in INTEGER_COLUMNS
        }
        actual_values = summarize_bin(unique_fasta_lengths(fasta_path), tiara_labels)
        for column in INTEGER_COLUMNS:
            if numeric_values[column] != actual_values[column]:
                raise ValueError(
                    f"Inconsistent {column} in bin-domain assignment row "
                    f"{line_number}: expected {actual_values[column]}, observed "
                    f"{numeric_values[column]}."
                )

        expected_domain = select_domain(
            {
                "prok": actual_values["prok_bp"],
                "euk": actual_values["euk_bp"],
            },
            fallback_domain,
        )
        if row["assigned_domain"] != expected_domain:
            raise ValueError(
                f"Inconsistent assigned_domain in bin-domain assignment row "
                f"{line_number}: expected {expected_domain!r}, observed "
                f"{row['assigned_domain']!r}."
            )

        validated_rows.append(
            {
                **row,
                **numeric_values,
                "source_fasta": fasta_path,
            }
        )
    return validated_rows


def main() -> None:
    refined_prok_dir = Path(snakemake.input.prok_bins)
    source_dirs = [Path(path) for path in snakemake.input.source_bins]
    assignments_path = Path(snakemake.input.bin_domain_assignments)
    tiara_tsv = Path(snakemake.input.tiara_tsv)
    final_prok_dir = Path(snakemake.output.prok_bins)
    final_euk_dir = Path(snakemake.output.euk_bins)
    final_unknown_dir = Path(snakemake.output.unknown_bins)
    summary_path = Path(snakemake.output.summary)
    annotation_samples_path = Path(snakemake.output.annotation_samples)
    rename_map_path = Path(snakemake.output.rename_map)
    log_path = Path(snakemake.log[0])
    tools = [str(tool) for tool in snakemake.params.tools]
    euk_source_tool = str(snakemake.params.euk_source_tool)
    fallback_domain = str(snakemake.params.fallback_domain)
    sample_prefix = str(snakemake.wildcards.sample)
    workdir = Path.cwd()

    if euk_source_tool not in tools:
        raise ValueError(
            f"Eukaryotic source tool '{euk_source_tool}' is not selected: {tools}."
        )

    tiara_labels = read_tiara_labels(tiara_tsv, fallback_domain)
    source_assignments = validated_assignments(
        assignments_path,
        tools,
        source_dirs,
        sample_prefix,
        fallback_domain,
        tiara_labels,
        workdir,
    )
    euk_sources = [
        row
        for row in source_assignments
        if row["tool"] == euk_source_tool and row["assigned_domain"] == "euk"
    ]
    unknown_sources = [
        row
        for row in source_assignments
        if row["tool"] == euk_source_tool and row["assigned_domain"] == "unknown"
    ]

    for output_dir in (final_prok_dir, final_euk_dir, final_unknown_dir):
        if output_dir.exists():
            shutil.rmtree(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
    for output_path in (
        summary_path,
        annotation_samples_path,
        rename_map_path,
        log_path,
    ):
        output_path.parent.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, str | int]] = []
    annotation_rows: list[dict[str, str]] = []
    rename_rows: list[dict[str, str]] = []
    domain_counts = defaultdict(int)

    final_sources: dict[str, list[tuple[Path, dict[str, int] | None]]] = {
        "prok": [
            (fasta_path, None) for fasta_path in direct_child_fastas(refined_prok_dir)
        ],
        "euk": [
            (
                row["source_fasta"],
                {column: row[column] for column in INTEGER_COLUMNS},
            )
            for row in euk_sources
        ],
        "unknown": [
            (
                row["source_fasta"],
                {column: row[column] for column in INTEGER_COLUMNS},
            )
            for row in unknown_sources
        ],
    }

    domain_dirs = {
        "prok": final_prok_dir,
        "euk": final_euk_dir,
        "unknown": final_unknown_dir,
    }
    for assigned_domain in ("prok", "euk", "unknown"):
        target_dir = domain_dirs[assigned_domain]
        for fasta_path, prepared_statistics in final_sources[assigned_domain]:
            statistics = prepared_statistics or summarize_bin(
                unique_fasta_lengths(fasta_path), tiara_labels
            )
            domain_counts[assigned_domain] += 1
            sample_name = generated_mag_identifier(
                sample_prefix, assigned_domain, domain_counts[assigned_domain]
            )
            standardized_name = f"{sample_name}.fa"
            standardized_path = target_dir / standardized_name
            shutil.copyfile(fasta_path, standardized_path)

            original_path = serialize_provenance_path(fasta_path, workdir)
            standardized_relative_path = serialize_provenance_path(
                standardized_path, workdir
            )
            rows.append(
                {
                    "bin": fasta_path.name,
                    "original_bin": fasta_path.name,
                    "original_path": original_path,
                    "standardized_bin": standardized_name,
                    "standardized_path": standardized_relative_path,
                    "assigned_domain": assigned_domain,
                    **statistics,
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

    summary_columns = [
        "bin",
        "original_bin",
        "original_path",
        "standardized_bin",
        "standardized_path",
        "assigned_domain",
        *INTEGER_COLUMNS,
        "fallback_domain",
    ]
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=summary_columns, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)

    annotation_columns = ["sample", "path", "domain", "genome_type"]
    with annotation_samples_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=annotation_columns, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(annotation_rows)

    rename_columns = [
        "original_bin",
        "original_path",
        "standardized_bin",
        "standardized_path",
        "domain",
    ]
    with rename_map_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=rename_columns, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rename_rows)

    log_path.write_text(
        f"final_prok_bins={domain_counts['prok']}\n"
        f"final_euk_bins={domain_counts['euk']}\n"
        f"final_unknown_bins={domain_counts['unknown']}\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
