"""Snakemake local wrapper for SemiBin2 with stable downstream outputs."""

import csv
import shlex
import shutil
from pathlib import Path

from snakemake.shell import shell

RESERVED_LONG_OPTIONS = {
    "--abundance",
    "--abundances",
    "--biome",
    "--compression",
    "--depth-metabat2",
    "--environment",
    "--habitat",
    "--input-bam",
    "--input-fasta",
    "--model",
    "--no-self-supervised",
    "--no-semi-supervised",
    "--output",
    "--processes",
    "--self-supervised",
    "--semi-supervised",
    "--sequencing-type",
    "--threads",
    "--training-type",
}
RESERVED_SHORT_OPTIONS = {"-a", "-b", "-i", "-o", "-p", "-t"}
BINS_INFO_COLUMNS = ["filename", "nbps", "nr_contigs", "N50", "L50"]


def _reserved_extra_option(token):
    if token == "--":
        return token
    option = token.split("=", 1)[0]
    if option in RESERVED_LONG_OPTIONS:
        return option
    if option.startswith("--") and any(
        reserved.startswith(option) for reserved in RESERVED_LONG_OPTIONS
    ):
        return option
    for short_option in RESERVED_SHORT_OPTIONS:
        if token == short_option or token.startswith(short_option) and len(token) > 2:
            return short_option
    return None


def _validated_extra(extra):
    try:
        extra_args = shlex.split(str(extra), posix=True)
    except ValueError as error:
        raise ValueError(
            "config.binning.semibin2.extra must be a valid shell-style option "
            f"string: {error}"
        ) from error

    reserved = sorted(
        {
            reserved
            for token in extra_args
            if (reserved := _reserved_extra_option(token)) is not None
        }
    )
    if reserved:
        raise ValueError(
            "config.binning.semibin2.extra must not set wrapper-managed "
            f"argument(s): {', '.join(reserved)}. Configure BAMs, input FASTA, "
            "output, coverage input, model mode/environment, sequencing type, "
            "compression, and threads through their declared MeeW inputs and "
            "settings."
        )
    return extra_args


def _canonical_bin_reference(filename, canonical_bins):
    source = Path(filename)
    bin_name = source.name
    if not bin_name:
        raise ValueError("SemiBin2 bins_info.tsv contains an empty filename value.")
    canonical_bin = canonical_bins / bin_name
    if not canonical_bin.is_file():
        raise FileNotFoundError(
            "SemiBin2 bins metadata references a bin that is absent from the "
            f"canonical output directory: {canonical_bin}"
        )
    if source.is_absolute():
        return str(canonical_bin.resolve())
    return str(canonical_bin)


def _normalize_bins_info(canonical_info, canonical_bins):
    if not canonical_info.is_file():
        if any(canonical_bins.iterdir()):
            raise FileNotFoundError(
                "SemiBin2 produced bins in output_bins/ but did not produce "
                "bins_info.tsv or recluster_bins_info.tsv."
            )
        canonical_info.write_text("\t".join(BINS_INFO_COLUMNS) + "\n", encoding="utf-8")
        return

    with canonical_info.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "filename" not in reader.fieldnames:
            raise ValueError(
                "SemiBin2 bins metadata must contain a 'filename' column: "
                f"{canonical_info}"
            )
        fieldnames = reader.fieldnames
        rows = list(reader)

    for row in rows:
        row["filename"] = _canonical_bin_reference(
            row.get("filename", ""), canonical_bins
        )

    temporary_info = canonical_info.with_name(f".{canonical_info.name}.tmp")
    try:
        with temporary_info.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)
        temporary_info.replace(canonical_info)
    finally:
        temporary_info.unlink(missing_ok=True)


output_dir = Path(str(snakemake.params.output_dir))
if output_dir.exists():
    shutil.rmtree(output_dir)
output_dir.parent.mkdir(parents=True, exist_ok=True)

if snakemake.log:
    Path(str(snakemake.log[0])).parent.mkdir(parents=True, exist_ok=True)
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

raw_bams = snakemake.input.bams
bams = [str(raw_bams)] if isinstance(raw_bams, str) else list(raw_bams)
if not bams:
    raise ValueError("SemiBin2 requires at least one resolved BAM input.")

environment = snakemake.params.environment
if len(bams) > 1 or environment is None:
    model_mode = ["--self-supervised"]
else:
    model_mode = ["--environment", str(environment)]
extra_args = _validated_extra(snakemake.params.extra)

shell(
    "SemiBin2 single_easy_bin"
    " -i {snakemake.input.contigs:q}"
    " -b {bams:q}"
    " -o {snakemake.params.output_dir:q}"
    " --threads {snakemake.threads}"
    " {extra_args:q}"
    " --sequencing-type {snakemake.params.sequencing_type:q}"
    " --compression none"
    " {model_mode:q}"
    " {log}"
)

canonical_bins = Path(str(snakemake.output.bins))
reclustered_bins = output_dir / "output_recluster_bins"
if reclustered_bins.is_dir():
    if canonical_bins.exists():
        shutil.rmtree(canonical_bins)
    reclustered_bins.replace(canonical_bins)
if not canonical_bins.is_dir():
    raise FileNotFoundError(
        "SemiBin2 did not produce output_bins/ or output_recluster_bins/."
    )

canonical_info = Path(str(snakemake.output.bins_info))
reclustered_info = output_dir / "recluster_bins_info.tsv"
if reclustered_info.is_file():
    if canonical_info.exists():
        canonical_info.unlink()
    reclustered_info.replace(canonical_info)
_normalize_bins_info(canonical_info, canonical_bins)
