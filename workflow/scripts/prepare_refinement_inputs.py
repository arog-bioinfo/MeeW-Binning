"""Prepare sample-level DAS Tool mappings and bin-domain metadata."""

import csv
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
    read_tiara_labels,
    select_domain,
    serialize_provenance_path,
    unique_fasta_lengths,
)


def main() -> None:
    bins_dirs = [Path(path) for path in snakemake.input.bins_dirs]
    assembly = Path(snakemake.input.assembly)
    tiara_tsv = Path(snakemake.input.tiara_tsv)
    mapping_paths = [Path(path) for path in snakemake.output.mappings]
    assignments_path = Path(snakemake.output.assignments)
    log_path = Path(snakemake.log[0])
    tools = [str(tool) for tool in snakemake.params.tools]
    fallback_domain = str(snakemake.params.fallback_domain)
    sample = str(snakemake.wildcards.sample)

    if len(tools) != len(bins_dirs) or len(tools) != len(mapping_paths):
        raise ValueError(
            "Refinement tools, source bin directories, and DAS Tool mappings must "
            "have matching lengths."
        )
    if len(set(tools)) != len(tools):
        raise ValueError("Refinement tools must be unique and ordered.")
    for tool, mapping_path in zip(tools, mapping_paths, strict=True):
        expected_name = f"{tool}_prok_bins.tsv"
        if mapping_path.name != expected_name:
            raise ValueError(
                f"DAS Tool mapping for '{tool}' must be named {expected_name}, "
                f"observed {mapping_path.name}."
            )

    assembly_contigs = set(unique_fasta_lengths(assembly))
    tiara_labels = read_tiara_labels(tiara_tsv, fallback_domain)
    assignment_rows: list[dict[str, str | int]] = []
    mappings_by_tool: dict[str, list[tuple[str, str]]] = {}
    log_rows = []

    for tool, bins_dir in zip(tools, bins_dirs, strict=True):
        contig_owner: dict[str, str] = {}
        mapping_rows: list[tuple[str, str]] = []
        domain_counts = defaultdict(int)
        bin_fastas = direct_child_fastas(bins_dir)

        for index, fasta_path in enumerate(bin_fastas, start=1):
            lengths = unique_fasta_lengths(fasta_path)
            for contig_id in lengths:
                if contig_id not in assembly_contigs:
                    raise ValueError(
                        f"Binner '{tool}' bin {fasta_path.name} contains contig "
                        f"'{contig_id}' that is absent from assembly {assembly}."
                    )
                if contig_id in contig_owner:
                    raise ValueError(
                        f"Binner '{tool}' assigns contig '{contig_id}' to both "
                        f"{contig_owner[contig_id]} and {fasta_path.name}."
                    )
                contig_owner[contig_id] = fasta_path.name

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
            domain_counts[assigned_domain] += 1
            normalized_bin_id = f"{sample}__{tool}__bin_{index:06d}"
            if assigned_domain == "prok":
                mapping_rows.extend(
                    (contig_id, normalized_bin_id) for contig_id in sorted(lengths)
                )

            assignment_rows.append(
                {
                    "tool": tool,
                    "source_bin": fasta_path.name,
                    "source_path": serialize_provenance_path(fasta_path, Path.cwd()),
                    "normalized_bin_id": normalized_bin_id,
                    "assigned_domain": assigned_domain,
                    "prok_bp": bp_by_domain.get("prok", 0),
                    "euk_bp": bp_by_domain.get("euk", 0),
                    "labeled_bp": labeled_bp,
                    "unlabeled_bp": unlabeled_bp,
                    "total_bp": labeled_bp + unlabeled_bp,
                    "n_contigs": len(lengths),
                    "n_labeled_contigs": sum(
                        1 for contig_id in lengths if contig_id in tiara_labels
                    ),
                    "fallback_domain": fallback_domain,
                }
            )

        mappings_by_tool[tool] = mapping_rows
        log_rows.extend(
            [
                f"{tool}.source_bins={len(bin_fastas)}",
                f"{tool}.prok_bins={domain_counts['prok']}",
                f"{tool}.euk_bins={domain_counts['euk']}",
                f"{tool}.unknown_bins={domain_counts['unknown']}",
                f"{tool}.prok_contig_assignments={len(mapping_rows)}",
            ]
        )

    assignments_path.parent.mkdir(parents=True, exist_ok=True)
    for tool, mapping_path in zip(tools, mapping_paths, strict=True):
        mapping_path.parent.mkdir(parents=True, exist_ok=True)
        with mapping_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerows(mappings_by_tool[tool])

    with assignments_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=ASSIGNMENT_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(assignment_rows)

    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(
        "\n".join(log_rows) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
