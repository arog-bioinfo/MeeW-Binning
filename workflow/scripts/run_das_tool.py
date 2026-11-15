"""Run DAS Tool 1.1.7 and normalize selected and valid-empty outcomes."""

import csv
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

_SCRIPTS_DIR = (
    Path(snakemake.scriptdir).resolve()
    if "snakemake" in globals() and hasattr(snakemake, "scriptdir")
    else Path(__file__).resolve().parent
)
sys.path.insert(0, str(_SCRIPTS_DIR))

from meew_binning_common import ASSIGNMENT_COLUMNS  # noqa: E402

SUMMARY_COLUMNS = [
    "bin",
    "bin_set",
    "unique_SCGs",
    "redundant_SCGs",
    "SCG_set",
    "size",
    "contigs",
    "N50",
    "bin_score",
    "SCG_completeness",
    "SCG_redundancy",
]
VALID_EMPTY_MESSAGES = (
    re.compile(r"No bins with bin-score >.* found\.", re.IGNORECASE),
    re.compile(r"No single copy genes predicted", re.IGNORECASE),
)
DOMAIN_COUNT_COLUMNS = {
    "prok": "prokaryotic_bins",
    "euk": "eukaryotic_bins",
    "unknown": "unknown_bins",
}


def write_header_only(path: Path, columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)


def mapping_counts(path: Path) -> tuple[int, set[str]]:
    contigs = 0
    bins: set[str] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, row in enumerate(reader, start=1):
            if not row:
                continue
            if len(row) != 2 or not row[0] or not row[1]:
                raise ValueError(
                    f"Invalid DAS Tool contig-to-bin row {line_number} in {path}."
                )
            contigs += 1
            bins.add(row[1])
    return contigs, bins


def assignment_counts(
    path: Path, labels: list[str]
) -> tuple[dict[str, dict[str, int]], dict[str, set[str]]]:
    """Count assignments and collect expected non-empty prokaryotic bin IDs."""
    counts = {
        label: {
            "source_bins": 0,
            "prokaryotic_bins": 0,
            "eukaryotic_bins": 0,
            "unknown_bins": 0,
        }
        for label in labels
    }
    prokaryotic_ids = {label: set() for label in labels}
    observed_ids: set[str] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ASSIGNMENT_COLUMNS:
            raise ValueError(
                f"Invalid bin-domain assignment columns in {path}: "
                f"expected {ASSIGNMENT_COLUMNS}, observed {reader.fieldnames}."
            )
        for line_number, row in enumerate(reader, start=2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(
                    f"Malformed bin-domain assignment row {line_number} in {path}."
                )
            tool = row["tool"]
            domain = row["assigned_domain"]
            if tool not in counts:
                raise ValueError(
                    f"Unknown tool '{tool}' in bin-domain assignment row "
                    f"{line_number} in {path}."
                )
            if domain not in DOMAIN_COUNT_COLUMNS:
                raise ValueError(
                    f"Invalid assigned_domain '{domain}' in row {line_number} of {path}."
                )
            normalized_bin_id = row["normalized_bin_id"]
            if not normalized_bin_id or normalized_bin_id in observed_ids:
                raise ValueError(
                    f"Missing or duplicate normalized_bin_id {normalized_bin_id!r} "
                    f"in row {line_number} of {path}."
                )
            observed_ids.add(normalized_bin_id)
            try:
                n_contigs = int(row["n_contigs"])
            except ValueError as error:
                raise ValueError(
                    f"Invalid n_contigs {row['n_contigs']!r} in row {line_number} "
                    f"of {path}."
                ) from error
            if n_contigs < 0 or str(n_contigs) != row["n_contigs"]:
                raise ValueError(
                    f"Invalid n_contigs {row['n_contigs']!r} in row {line_number} "
                    f"of {path}."
                )
            counts[tool]["source_bins"] += 1
            counts[tool][DOMAIN_COUNT_COLUMNS[domain]] += 1
            if domain == "prok" and n_contigs > 0:
                prokaryotic_ids[tool].add(normalized_bin_id)
    return counts, prokaryotic_ids


def copy_or_header(source: Path, destination: Path, columns: list[str]) -> None:
    if source.is_file() and source.stat().st_size > 0:
        shutil.copyfile(source, destination)
    else:
        write_header_only(destination, columns)


def write_status(
    path: Path,
    *,
    status: str,
    reason: str,
    version: str,
    search_engine: str,
    score_threshold: float,
    threads: int,
    selected_bins: int,
) -> None:
    columns = [
        "status",
        "reason",
        "das_tool_version",
        "search_engine",
        "score_threshold",
        "threads",
        "selected_bins",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerow(
            {
                "status": status,
                "reason": reason,
                "das_tool_version": version,
                "search_engine": search_engine,
                "score_threshold": score_threshold,
                "threads": threads,
                "selected_bins": selected_bins,
            }
        )


def make_empty_outputs(
    *,
    bins_dir: Path,
    summary: Path,
    scores: Path,
    contig2bin: Path,
    tool_log: Path,
    tool_log_text: str,
) -> None:
    if bins_dir.exists():
        shutil.rmtree(bins_dir)
    bins_dir.mkdir(parents=True, exist_ok=True)
    write_header_only(summary, SUMMARY_COLUMNS)
    write_header_only(scores, SUMMARY_COLUMNS)
    contig2bin.write_text("", encoding="utf-8")
    tool_log.write_text(tool_log_text.rstrip() + "\n", encoding="utf-8")


def main() -> None:
    assembly = Path(snakemake.input.assembly)
    mapping_paths = [Path(path) for path in snakemake.input.contig2bin]
    bin_domain_assignments = Path(snakemake.input.bin_domain_assignments)
    labels = [str(label) for label in snakemake.params.labels]
    expected_version = str(snakemake.params.expected_version)
    search_engine = str(snakemake.params.search_engine)
    score_threshold = float(snakemake.params.score_threshold)
    threads = int(snakemake.threads)

    bins_dir = Path(snakemake.output.bins)
    summary = Path(snakemake.output.summary)
    scores = Path(snakemake.output.scores)
    contig2bin = Path(snakemake.output.contig2bin)
    provenance = Path(snakemake.output.provenance)
    run_status = Path(snakemake.output.run_status)
    tool_log = Path(snakemake.output.tool_log)
    execution_log = Path(snakemake.log[0])

    if len(mapping_paths) != len(labels):
        raise ValueError(
            "DAS Tool labels and contig-to-bin inputs must have matching lengths."
        )
    for label, mapping_path in zip(labels, mapping_paths, strict=True):
        expected_name = f"{label}_prok_bins.tsv"
        if mapping_path.name != expected_name:
            raise ValueError(
                f"DAS Tool mapping for '{label}' must be named {expected_name}, "
                f"observed {mapping_path.name}."
            )
    if search_engine != "diamond":
        raise ValueError("MeeW-Binning requires DAS Tool search_engine='diamond'.")
    if assembly.suffix.lower() == ".gz":
        raise ValueError("DAS Tool 1.1.7 requires an uncompressed assembly FASTA.")

    output_parent = summary.parent
    output_parent.mkdir(parents=True, exist_ok=True)
    execution_log.parent.mkdir(parents=True, exist_ok=True)
    for output_path in (summary, scores, contig2bin, provenance, run_status, tool_log):
        output_path.parent.mkdir(parents=True, exist_ok=True)

    version_result = subprocess.run(
        ["DAS_Tool", "--version"],
        check=True,
        capture_output=True,
        text=True,
    )
    version_text = (version_result.stdout + version_result.stderr).strip()
    if expected_version not in version_text:
        raise RuntimeError(
            f"Expected DAS Tool {expected_version}, observed: {version_text!r}."
        )

    total_contigs = 0
    provenance_rows = []
    assignment_bin_counts, expected_prokaryotic_ids = assignment_counts(
        bin_domain_assignments, labels
    )
    for label, mapping_path in zip(labels, mapping_paths, strict=True):
        contig_count, mapped_bin_ids = mapping_counts(mapping_path)
        if mapped_bin_ids != expected_prokaryotic_ids[label]:
            missing_ids = sorted(expected_prokaryotic_ids[label] - mapped_bin_ids)
            unexpected_ids = sorted(mapped_bin_ids - expected_prokaryotic_ids[label])
            raise ValueError(
                f"DAS Tool mapping IDs for '{label}' do not match consolidated "
                f"prokaryotic assignments; missing={missing_ids}, "
                f"unexpected={unexpected_ids}."
            )
        total_contigs += contig_count
        provenance_rows.append(
            {
                "label": label,
                "contig2bin": str(mapping_path),
                "bin_domain_assignments": str(bin_domain_assignments),
                **assignment_bin_counts[label],
                "prokaryotic_contigs": contig_count,
                "mapped_prokaryotic_bins": len(mapped_bin_ids),
            }
        )

    provenance_columns = [
        "label",
        "contig2bin",
        "bin_domain_assignments",
        "source_bins",
        "prokaryotic_bins",
        "eukaryotic_bins",
        "unknown_bins",
        "prokaryotic_contigs",
        "mapped_prokaryotic_bins",
    ]
    with provenance.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=provenance_columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(provenance_rows)

    if total_contigs == 0:
        reason = "No prokaryotic DAS Tool input contigs were available."
        make_empty_outputs(
            bins_dir=bins_dir,
            summary=summary,
            scores=scores,
            contig2bin=contig2bin,
            tool_log=tool_log,
            tool_log_text=(
                f"DAS Tool {expected_version}\n"
                "MeeW-Binning status: no_prokaryotic_inputs\n"
                f"{reason}\nThe DAS Tool refinement command was not invoked."
            ),
        )
        execution_log.write_text(reason + "\n", encoding="utf-8")
        write_status(
            run_status,
            status="empty",
            reason="no_prokaryotic_inputs",
            version=expected_version,
            search_engine=search_engine,
            score_threshold=score_threshold,
            threads=threads,
            selected_bins=0,
        )
        return

    if bins_dir.exists():
        shutil.rmtree(bins_dir)
    bins_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(
        prefix=".das_tool.", dir=output_parent
    ) as temporary_directory:
        prefix = Path(temporary_directory) / "run"
        command = [
            "DAS_Tool",
            "--bins",
            ",".join(str(path) for path in mapping_paths),
            "--labels",
            ",".join(labels),
            "--contigs",
            str(assembly),
            "--outputbasename",
            str(prefix),
            "--search_engine",
            search_engine,
            "--write_bin_evals",
            "--write_bins",
            "--score_threshold",
            str(score_threshold),
            "--threads",
            str(threads),
            "--debug",
        ]
        with execution_log.open("w", encoding="utf-8") as log_handle:
            log_handle.write("Command: " + " ".join(command) + "\n")
            log_handle.flush()
            result = subprocess.run(
                command,
                check=False,
                stdout=log_handle,
                stderr=subprocess.STDOUT,
                text=True,
            )

        raw_tool_log = Path(f"{prefix}_DASTool.log")
        diagnostic_text = execution_log.read_text(encoding="utf-8")
        if raw_tool_log.is_file():
            diagnostic_text += "\n" + raw_tool_log.read_text(
                encoding="utf-8", errors="replace"
            )

        if result.returncode != 0:
            empty_reason = next(
                (
                    pattern.pattern
                    for pattern in VALID_EMPTY_MESSAGES
                    if pattern.search(diagnostic_text)
                ),
                None,
            )
            if empty_reason is None:
                raise subprocess.CalledProcessError(result.returncode, command)

            raw_scores = Path(f"{prefix}_allBins.eval")
            make_empty_outputs(
                bins_dir=bins_dir,
                summary=summary,
                scores=scores,
                contig2bin=contig2bin,
                tool_log=tool_log,
                tool_log_text=(
                    raw_tool_log.read_text(encoding="utf-8", errors="replace")
                    if raw_tool_log.is_file()
                    else diagnostic_text
                ),
            )
            copy_or_header(raw_scores, scores, SUMMARY_COLUMNS)
            write_status(
                run_status,
                status="empty",
                reason="no_das_tool_selected_bins",
                version=expected_version,
                search_engine=search_engine,
                score_threshold=score_threshold,
                threads=threads,
                selected_bins=0,
            )
            return

        raw_summary = Path(f"{prefix}_DASTool_summary.tsv")
        raw_scores = Path(f"{prefix}_allBins.eval")
        raw_contig2bin = Path(f"{prefix}_DASTool_contig2bin.tsv")
        raw_bins = Path(f"{prefix}_DASTool_bins")
        required_outputs = [
            raw_summary,
            raw_scores,
            raw_contig2bin,
            raw_bins,
            raw_tool_log,
        ]
        missing = [str(path) for path in required_outputs if not path.exists()]
        if missing:
            raise RuntimeError(
                "DAS Tool completed without required output(s): " + ", ".join(missing)
            )

        selected_fastas = sorted(
            path
            for path in raw_bins.iterdir()
            if path.is_file() and path.suffix.lower() in {".fa", ".fna", ".fasta"}
        )
        if not selected_fastas:
            raise RuntimeError(
                "DAS Tool completed successfully but selected no FASTA bins."
            )

        shutil.copyfile(raw_summary, summary)
        shutil.copyfile(raw_scores, scores)
        shutil.copyfile(raw_contig2bin, contig2bin)
        shutil.copyfile(raw_tool_log, tool_log)
        for fasta_path in selected_fastas:
            shutil.copyfile(fasta_path, bins_dir / fasta_path.name)

        write_status(
            run_status,
            status="selected",
            reason="selected_bins",
            version=expected_version,
            search_engine=search_engine,
            score_threshold=score_threshold,
            threads=threads,
            selected_bins=len(selected_fastas),
        )


if __name__ == "__main__":
    main()
