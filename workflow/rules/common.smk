import csv
import importlib.util
import os
from pathlib import Path

import pandas as pd
from snakemake.exceptions import WorkflowError
from snakemake.utils import validate

SAMPLES_SCHEMA = workflow.source_path("../schemas/samples.schema.yaml")
CHECKM1_HELPER = workflow.source_path("../scripts/common/checkm1.py")
CHECKM1_RUNNER = workflow.source_path("../scripts/run_comebin_with_checkm.py")
# Compatibility aliases retained for Main configurations/tests that imported the
# original COMEBin-specific names before the runtime became tool-aware.
COMEBIN_CHECKM1_HELPER = CHECKM1_HELPER
COMEBIN_CHECKM1_RUNNER = CHECKM1_RUNNER


def _load_dag_helper(module_name, source):
    spec = importlib.util.spec_from_file_location(module_name, source)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load Binning helper from {source}")
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    return helper


_checkm1_helpers = _load_dag_helper("_meew_binning_checkm1_paths", Path(CHECKM1_HELPER))


def expand_path(value):
    """Expand environment variables and ``~`` while preserving relative paths."""
    if value is None or str(value).strip() == "":
        return ""
    return os.path.expanduser(os.path.expandvars(str(value)))


# --- Main config validation -------------------------------------------------

validate(config, schema="../schemas/config.schema.yaml")

RESOURCES_ROOT = expand_path(config.get("resources_root", ""))
CHECKM1_CONFIG = config.get("checkm1", {})
CHECKM1_DATABASE_PATH = _checkm1_helpers.resolve_checkm1_database(
    CHECKM1_CONFIG.get("database_path", ""),
    RESOURCES_ROOT,
)
CHECKM1_DOWNLOAD_DATABASE = CHECKM1_CONFIG.get("download_database", True)

metabinner_config = config.get("binning", {}).get("metabinner", {})
METABINNER_CONFIG = metabinner_config
metabinner_config["install_path"] = expand_path(
    metabinner_config.get("install_path", "")
)

coverm_config = config.get("coverm", {})
for field in ["bam", "reads_1", "reads_2", "single", "interleaved"]:
    if field not in coverm_config:
        continue
    if isinstance(coverm_config[field], dict):
        coverm_config[field] = {
            sample: expand_path(value) for sample, value in coverm_config[field].items()
        }
    else:
        coverm_config[field] = expand_path(coverm_config[field])

# --- Sample sheet loading ----------------------------------------------------

_samples_cache = None


def _configured_sample_sheet_path():
    value = config.get("sample_sheet")
    if value is None or str(value).strip() == "":
        raise WorkflowError(
            "Binning sample metadata requires a non-empty `sample_sheet` "
            "configuration value."
        )
    path = expand_path(value)
    config["sample_sheet"] = path
    return path


def load_samples():
    global _samples_cache
    if _samples_cache is None:
        sample_sheet = Path(_configured_sample_sheet_path())
        if not sample_sheet.exists():
            raise WorkflowError(f"Sample sheet does not exist: {sample_sheet}")

        try:
            with sample_sheet.open(newline="", encoding="utf-8") as handle:
                reader = csv.DictReader(handle, delimiter="\t", strict=True)
                fieldnames = [field.strip() for field in (reader.fieldnames or [])]
                if not fieldnames or any(not field for field in fieldnames):
                    raise WorkflowError(
                        f"Sample sheet must have a non-empty TSV header: {sample_sheet}"
                    )

                duplicate_headers = sorted(
                    {field for field in fieldnames if fieldnames.count(field) > 1}
                )
                if duplicate_headers:
                    raise WorkflowError(
                        "Sample sheet contains duplicate header(s): "
                        + ", ".join(duplicate_headers)
                    )
                reader.fieldnames = fieldnames

                rows = []
                for row in reader:
                    if None in row:
                        raise WorkflowError(
                            f"Sample sheet row {reader.line_num} has more fields than the header"
                        )
                    if any(value is None for value in row.values()):
                        raise WorkflowError(
                            f"Sample sheet row {reader.line_num} has fewer fields than the header"
                        )

                    raw_row = {
                        key: value if key == "sample" else value.strip()
                        for key, value in row.items()
                    }
                    validate(raw_row, schema=SAMPLES_SCHEMA)

                    # Blank optional cells participate in schema validation and
                    # are expanded only afterward when they contain a path.
                    normalized = dict(raw_row)
                    for column in fieldnames:
                        if column != "sample" and normalized[column]:
                            normalized[column] = expand_path(normalized[column])
                    rows.append(normalized)
        except csv.Error as error:
            raise WorkflowError(
                f"Malformed TSV sample sheet {sample_sheet}: {error}"
            ) from error

        if not rows:
            raise WorkflowError(
                f"Sample sheet must contain at least one row: {sample_sheet}"
            )

        sample_ids = [row["sample"] for row in rows]
        duplicate_samples = sorted(
            {sample for sample in sample_ids if sample_ids.count(sample) > 1}
        )
        if duplicate_samples:
            raise WorkflowError(
                "Sample sheet contains duplicate sample identifier(s): "
                + ", ".join(duplicate_samples)
            )

        _samples_cache = (
            pd.DataFrame(rows, columns=fieldnames)
            .set_index("sample", drop=False)
            .sort_index()
        )
        _validate_config_dependent_sample_fields(_samples_cache)
    return _samples_cache


def _sample_value(samples_df, sample, field):
    if field not in samples_df.columns:
        return ""
    return str(samples_df.loc[sample, field])


def _coverm_value(samples_df, sample, field):
    if field not in COVERM_CONFIG:
        return _sample_value(samples_df, sample, field)
    configured = COVERM_CONFIG[field]
    if isinstance(configured, dict):
        return str(configured.get(sample, ""))
    return str(configured)


def _validate_config_dependent_sample_fields(samples_df):
    binner_fields = {
        "comebin": ["bam"],
        "metabinner": ["bam"],
        "semibin2": ["bam"],
        "metabat2": ["bam"],
    }
    for sample in samples_df.index:
        for selected_binner in SELECTED_BINNING_TOOLS:
            for field in binner_fields[selected_binner]:
                if not _sample_value(samples_df, sample, field).strip():
                    raise WorkflowError(
                        f"Sample '{sample}' requires non-empty '{field}' for binner "
                        f"'{selected_binner}'."
                    )

        coverm_fields = {
            "bam": ["bam"],
            "coupled": ["reads_1", "reads_2"],
            "single": ["single"],
            "interleaved": ["interleaved"],
        }
        for field in coverm_fields[COVERM_INPUT_MODE]:
            if not _coverm_value(samples_df, sample, field).strip():
                raise WorkflowError(
                    f"Sample '{sample}' requires non-empty '{field}' for "
                    f"CoverM input_mode '{COVERM_INPUT_MODE}'."
                )


# --- Constants and allowed tool sets -----------------------------------------

RESULTS_ROOT = "results/{sample}"
LOGS_ROOT = "<results>/logs"
RESULTS_DIR = "<results>"

ALLOWED_BINNING_TOOLS = {"comebin", "metabinner", "semibin2", "metabat2"}


# --- Tool-specific configuration ---------------------------------------------

SELECTED_BINNING_TOOLS = config["binning"]["tools"]
if not SELECTED_BINNING_TOOLS:
    raise WorkflowError("config.binning.tools must contain at least one binning tool.")
if not set(SELECTED_BINNING_TOOLS).issubset(ALLOWED_BINNING_TOOLS):
    raise WorkflowError(
        f"Unsupported binning tool configured: {set(SELECTED_BINNING_TOOLS) - ALLOWED_BINNING_TOOLS}"
    )

REFINEMENT_CONFIG = config["refinement"]
REFINEMENT_ENABLED = REFINEMENT_CONFIG["enabled"]
REFINEMENT_METHOD = REFINEMENT_CONFIG["method"]
EUK_SOURCE_TOOL = REFINEMENT_CONFIG.get("euk_source_tool", "")
DASTOOL_CONFIG = REFINEMENT_CONFIG.get("das_tool", {})

if REFINEMENT_ENABLED:
    if REFINEMENT_METHOD != "das_tool":
        raise WorkflowError(
            "config.refinement.method must be 'das_tool' when refinement is enabled."
        )
    if len(SELECTED_BINNING_TOOLS) < 2:
        raise WorkflowError(
            "DAS Tool refinement requires at least two distinct binning.tools."
        )
    if EUK_SOURCE_TOOL not in SELECTED_BINNING_TOOLS:
        raise WorkflowError(
            "config.refinement.euk_source_tool must be one of the selected "
            "config.binning.tools."
        )
else:
    if REFINEMENT_METHOD != "none":
        raise WorkflowError(
            "config.refinement.method must be 'none' when refinement is disabled."
        )
    if len(SELECTED_BINNING_TOOLS) != 1:
        raise WorkflowError(
            "Refinement-disabled mode requires exactly one binning tool."
        )

DASTOOL_THREADS = DASTOOL_CONFIG.get("threads", 8)
DASTOOL_SCORE_THRESHOLD = DASTOOL_CONFIG.get("score_threshold", 0.5)

COVERM_CONFIG = config.get("coverm", {})
COVERM_DOMAINS = COVERM_CONFIG.get("domains", ["prok", "euk"])
COVERM_INPUT_MODE = COVERM_CONFIG.get("input_mode", "bam")
COVERM_MAPPER = COVERM_CONFIG.get("mapper", "minimap2-sr")
COVERM_METHODS = COVERM_CONFIG.get(
    "methods", ["relative_abundance", "mean", "covered_fraction", "length", "count"]
)
COVERM_THREADS = COVERM_CONFIG.get("threads", 4)
COVERM_EXTRA = COVERM_CONFIG.get("extra", "")
COVERM_GENOME_FASTA_EXTENSION = COVERM_CONFIG.get("genome_fasta_extension", ".fa")


# --- Configuration validation ------------------------------------------------

FALLBACK_DOMAIN = config["domain_separation"].get("fallback_domain", "unknown")
if FALLBACK_DOMAIN not in {"prok", "euk", "unknown"}:
    raise ValueError(
        "config.domain_separation.fallback_domain must be one of: prok, euk, unknown"
    )

if COVERM_INPUT_MODE not in {"bam", "coupled", "single", "interleaved"}:
    raise ValueError(
        "config.coverm.input_mode must be one of: bam, coupled, single, interleaved"
    )


# --- Input helper functions --------------------------------------------------


def sample_input(wildcards, field):
    return str(load_samples().loc[wildcards.sample, field])


def comebin_bam_dir(bam_input):
    bam_path = Path(str(bam_input))
    if bam_path.suffix == ".bam":
        return str(bam_path.parent)
    return str(bam_path)


def resolve_bams(bam_input, sample):
    bam_path = Path(str(bam_input))
    if bam_path.is_dir():
        bam_files = sorted(path for path in bam_path.glob("*.bam") if path.is_file())
        if not bam_files:
            raise WorkflowError(
                f"Sample '{sample}' BAM directory contains no direct-child *.bam "
                f"files: {bam_path}"
            )
        return [str(path) for path in bam_files]

    if bam_path.suffix.lower() == ".bam":
        return [str(bam_path)]

    if bam_path.exists():
        raise WorkflowError(
            f"Sample '{sample}' bam value must be a .bam file or a directory "
            f"containing direct-child *.bam files: {bam_path}"
        )
    raise WorkflowError(
        f"Sample '{sample}' bam value does not exist and is not a .bam file path: "
        f"{bam_path}"
    )


def bam_inputs(wildcards):
    return resolve_bams(sample_input(wildcards, "bam"), wildcards.sample)


def coverm_sample_input(wildcards, field):
    if field in COVERM_CONFIG:
        value = COVERM_CONFIG[field]
        if isinstance(value, dict):
            return str(value[wildcards.sample])
        return str(value)
    return sample_input(wildcards, field)


def coverm_inputs(wildcards):
    inputs = {
        "bins": _format_rule_output(
            getattr(rules.separate_bins_by_domain.output, f"{wildcards.domain}_bins"),
            sample=wildcards.sample,
        )
    }

    if COVERM_INPUT_MODE == "bam":
        inputs["bam"] = coverm_sample_input(wildcards, "bam")
    elif COVERM_INPUT_MODE == "coupled":
        inputs["reads_1"] = coverm_sample_input(wildcards, "reads_1")
        inputs["reads_2"] = coverm_sample_input(wildcards, "reads_2")
    elif COVERM_INPUT_MODE == "single":
        inputs["single"] = coverm_sample_input(wildcards, "single")
    elif COVERM_INPUT_MODE == "interleaved":
        inputs["interleaved"] = coverm_sample_input(wildcards, "interleaved")

    return inputs


# --- Binning helper functions ------------------------------------------------


def _format_rule_output(template, **wildcards):
    return str(template).format(**wildcards)


def validate_selected_binner(tool):
    if tool not in ALLOWED_BINNING_TOOLS:
        raise WorkflowError(f"Unsupported binning tool requested: {tool}")
    if tool not in SELECTED_BINNING_TOOLS:
        raise WorkflowError(
            f"Binning tool '{tool}' is supported but not selected in "
            "config.binning.tools."
        )


def binner_bins_dir(sample, tool):
    validate_selected_binner(tool)
    if tool == "comebin":
        return _format_rule_output(rules.comebin.output.bins, sample=sample)
    if tool == "metabinner":
        return _format_rule_output(rules.metabinner_bins.output.bins, sample=sample)
    if tool == "metabat2":
        return _format_rule_output(rules.metabat2.output.bins, sample=sample)
    if tool == "semibin2":
        return _format_rule_output(rules.semibin2.output.bins, sample=sample)
    raise AssertionError(f"Unhandled selected binner: {tool}")


def selected_binner_bins_dir(sample):
    return binner_bins_dir(sample, SELECTED_BINNING_TOOLS[0])


def binner_tracked_output(sample, tool):
    validate_selected_binner(tool)
    if tool == "comebin":
        return _format_rule_output(rules.comebin.output.result_tsv, sample=sample)
    if tool == "metabinner":
        return _format_rule_output(rules.metabinner_bins.output.bins, sample=sample)
    if tool == "metabat2":
        return _format_rule_output(rules.metabat2.output.bins, sample=sample)
    if tool == "semibin2":
        return _format_rule_output(rules.semibin2.output.bins, sample=sample)
    raise AssertionError(f"Unhandled selected binner: {tool}")


def selected_binner_tracked_output(sample):
    return binner_tracked_output(sample, SELECTED_BINNING_TOOLS[0])


def refinement_binner_bins_inputs(wildcards):
    return [binner_bins_dir(wildcards.sample, tool) for tool in SELECTED_BINNING_TOOLS]


# --- Workflow targets ---------------------------------------------------------


def final_targets(wildcards):
    current_samples = load_samples()
    targets = []
    for sample in current_samples.index:
        targets.append(
            _format_rule_output(
                rules.separate_bins_by_domain.output.summary, sample=sample
            )
        )

        targets.extend(
            _format_rule_output(
                rules.coverm_genome_coverage.output[0],
                sample=sample,
                domain=domain,
            )
            for domain in COVERM_DOMAINS
        )

    return targets


ALL_TARGETS = final_targets
