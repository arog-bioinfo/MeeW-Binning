"""Snakemake local wrapper for CoverM."""

from pathlib import Path

from snakemake.shell import shell


def _as_list(value):
    if isinstance(value, (list, tuple)):
        return list(value)
    return [value]


log = snakemake.log_fmt_shell(stdout=True, stderr=True)
output_file = Path(snakemake.output[0])
output_file.parent.mkdir(parents=True, exist_ok=True)
bins_dir = Path(snakemake.input.bins)
extension = str(snakemake.params.extension)

has_matching_bin = any(
    path.is_file() and path.name.endswith(extension) for path in bins_dir.iterdir()
)
if not has_matching_bin:
    output_file.write_text("", encoding="utf-8")
    log_file = Path(snakemake.log[0])
    log_file.parent.mkdir(parents=True, exist_ok=True)
    log_file.write_text(
        "CoverM skipped: no direct-child regular files matching genome FASTA "
        f"extension {extension!r} in {bins_dir}.\n",
        encoding="utf-8",
    )
else:
    input_mode = snakemake.params.input_mode
    methods = _as_list(snakemake.params.methods)

    if input_mode == "bam":
        mapping_option = "--bam-files"
        mapping_values = _as_list(snakemake.input.bam)
    elif input_mode == "coupled":
        mapping_option = "--coupled"
        mapping_values = _as_list(snakemake.input.reads_1) + _as_list(
            snakemake.input.reads_2
        )
    elif input_mode == "single":
        mapping_option = "--single"
        mapping_values = _as_list(snakemake.input.single)
    elif input_mode == "interleaved":
        mapping_option = "--interleaved"
        mapping_values = _as_list(snakemake.input.interleaved)
    else:
        raise ValueError(f"Unsupported CoverM input_mode: {input_mode}")

    mapper = "" if input_mode == "bam" else snakemake.params.mapper or ""
    mapper_option = "--mapper" if mapper else ""

    shell(
        "coverm genome "
        " {mapping_option} {mapping_values:q} "
        " --genome-fasta-directory {snakemake.input.bins:q} "
        " --genome-fasta-extension {snakemake.params.extension:q} "
        " --methods {methods:q} "
        " {mapper_option} {mapper:q} "
        " --threads {snakemake.threads} "
        " --output-file {output_file:q} "
        " {snakemake.params.extra} "
        " {log}"
    )
