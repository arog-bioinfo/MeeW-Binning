import json
import subprocess
from pathlib import Path

import pytest
import yaml
from jsonschema import Draft202012Validator, ValidationError

ROOT = Path(__file__).resolve().parents[1]
SNAKEFILE = ROOT / "workflow" / "Snakefile"


def _module_config(sample_sheet_marker=None):
    config = yaml.safe_load(
        (ROOT / "config" / "config.yaml").read_text(encoding="utf-8")
    )
    config.pop("sample_sheet")
    if sample_sheet_marker is not None:
        config["sample_sheet"] = sample_sheet_marker
    return config


def _run_harness(tmp_path, harness, config):
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config), encoding="utf-8")
    return subprocess.run(
        [
            "snakemake",
            "--snakefile",
            str(harness),
            "--directory",
            str(tmp_path),
            "--configfile",
            str(config_path),
            "--replace-workflow-config",
            "--nolock",
            "--cores",
            "1",
            "--dry-run",
        ],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )


def test_module_can_import_scientific_rule_without_requesting_sample_metadata(tmp_path):
    contigs = tmp_path / "assembly.fa"
    contigs.write_text(">contig\nACGT\n", encoding="utf-8")
    harness = tmp_path / "sample-sheet-integration.smk"
    harness.write_text(
        f"""module binning:
    snakefile: {json.dumps(str(SNAKEFILE))}
    config: config
    pathvars:
        results="results/{{sample}}",
    skip_validation: True

use rule * from binning exclude all as binning_*

use rule tiara_predict from binning as binning_tiara_predict with:
    input:
        contigs="assembly.fa",

rule integration_target:
    default_target: True
    input:
        "results/sample-sheet/domain/tiara_contigs.tsv",
""",
        encoding="utf-8",
    )

    result = _run_harness(tmp_path, harness, _module_config())

    assert result.returncode == 0, result.stdout + result.stderr
    assert "binning_tiara_predict" in result.stdout
    assert "sample_sheet" not in result.stdout + result.stderr


@pytest.mark.parametrize("sample_sheet", [None, "", "   "])
def test_missing_or_blank_metadata_fails_only_when_loading_is_requested(
    tmp_path, sample_sheet
):
    harness = tmp_path / "request-standalone-metadata.smk"
    harness.write_text(
        f"""module binning:
    snakefile: {json.dumps(str(SNAKEFILE))}
    config: config
    pathvars:
        results="results/{{sample}}",
    skip_validation: True

use rule all from binning as requested_binning_all
""",
        encoding="utf-8",
    )

    result = _run_harness(tmp_path, harness, _module_config(sample_sheet))

    assert result.returncode != 0
    assert (
        "Binning sample metadata requires a non-empty `sample_sheet` configuration value."
        in result.stdout + result.stderr
    )


def test_standalone_schema_still_requires_nonempty_sample_sheet():
    schema = yaml.safe_load(
        (ROOT / "workflow" / "schemas" / "config.schema.yaml").read_text(
            encoding="utf-8"
        )
    )
    assert "sample_sheet" in schema["required"]
    assert schema["properties"]["sample_sheet"]["minLength"] == 1
    assert schema["properties"]["sample_sheet"]["pattern"] == r"\S"


def test_standalone_schema_rejects_removed_stage_sheets_setting():
    schema = yaml.safe_load(
        (ROOT / "workflow" / "schemas" / "config.schema.yaml").read_text(
            encoding="utf-8"
        )
    )
    config = yaml.safe_load((ROOT / "config" / "config.yaml").read_text())
    config["stage_sheets"] = {"enabled": False}

    with pytest.raises(ValidationError):
        Draft202012Validator(schema).validate(config)
