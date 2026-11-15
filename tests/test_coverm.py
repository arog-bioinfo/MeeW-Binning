import ast
import copy
import importlib
import json
import runpy
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml
from jsonschema import Draft202012Validator, ValidationError

ROOT = Path(__file__).resolve().parents[1]
SNAKEFILE = ROOT / "workflow" / "Snakefile"
SCRIPTS = ROOT / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

import finalize_refined_domains  # noqa: E402
import meew_binning_common  # noqa: E402
import prepare_refinement_inputs  # noqa: E402
import run_das_tool  # noqa: E402


def _schema():
    return yaml.safe_load(
        (ROOT / "workflow" / "schemas" / "config.schema.yaml").read_text(
            encoding="utf-8"
        )
    )


def test_coverm_domains_schema_and_shipped_configurations():
    schema = _schema()
    domains_schema = schema["properties"]["coverm"]["properties"]["domains"]
    assert domains_schema["default"] == ["prok", "euk"]
    assert domains_schema["minItems"] == 1
    assert domains_schema["uniqueItems"] is True
    assert domains_schema["items"]["enum"] == ["prok", "euk"]

    validator = Draft202012Validator(schema)
    shipped_configs = []
    for relative_path in (
        "config/config.yaml",
        ".test/config/config.yaml",
        ".test/config/config.single.yaml",
    ):
        config = yaml.safe_load((ROOT / relative_path).read_text(encoding="utf-8"))
        validator.validate(config)
        shipped_configs.append(config)

    assert shipped_configs[0]["coverm"]["domains"] == ["prok", "euk"]
    assert shipped_configs[1]["coverm"]["domains"] == ["prok", "euk"]
    assert all("enabled" not in config.get("coverm", {}) for config in shipped_configs)

    defaulted_config = copy.deepcopy(shipped_configs[0])
    defaulted_config["coverm"].pop("domains")
    validator.validate(defaulted_config)

    for invalid_domains in ([], ["prok", "prok"], ["unknown"]):
        invalid_config = copy.deepcopy(shipped_configs[0])
        invalid_config["coverm"]["domains"] = invalid_domains
        with pytest.raises(ValidationError):
            validator.validate(invalid_config)

    enabled_config = copy.deepcopy(shipped_configs[0])
    enabled_config["coverm"]["enabled"] = True
    with pytest.raises(ValidationError):
        validator.validate(enabled_config)


def _coverm_config(tmp_path, domains):
    sample = "coverm-sample"
    contigs = tmp_path / "assembly.fa"
    bam = tmp_path / "mapping.bam"
    contigs.write_text(">contig\nACGT\n", encoding="utf-8")
    bam.touch()
    sample_sheet = tmp_path / "samples.tsv"
    sample_sheet.write_text(
        f"sample\tcontigs\tbam\n{sample}\t{contigs}\t{bam}\n",
        encoding="utf-8",
    )
    coverm = {
        "input_mode": "bam",
        "mapper": "minimap2-sr",
        "methods": ["mean"],
        "threads": 1,
        "extra": "",
    }
    if domains is not None:
        coverm["domains"] = domains
    config = {
        "sample_sheet": str(sample_sheet),
        "binning": {
            "tools": ["comebin"],
            "comebin": {"extra": ""},
            "metabinner": {
                "install_path": "",
                "min_contig_length": 1000,
                "scale": "large",
            },
            "semibin2": {
                "environment": "global",
                "sequencing_type": "short_read",
                "extra": "",
            },
            "metabat2": {"min_contig_length": 1500, "seed": 1, "extra": ""},
        },
        "refinement": {"enabled": False, "method": "none"},
        "domain_separation": {
            "fallback_domain": "unknown",
            "tiara": {"extra": ""},
        },
        "coverm": coverm,
    }
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config), encoding="utf-8")
    return config_path, sample


@pytest.mark.parametrize(
    ("domains", "expected_domains"),
    [(["euk"], ["euk"]), (None, ["prok", "euk"])],
)
def test_final_coverm_targets_follow_configured_or_default_domains(
    tmp_path, domains, expected_domains
):
    config_path, sample = _coverm_config(tmp_path, domains)
    result = subprocess.run(
        [
            "snakemake",
            "--snakefile",
            str(SNAKEFILE),
            "--configfile",
            str(config_path),
            "--directory",
            str(tmp_path),
            "--cores",
            "1",
            "--nolock",
            "--dry-run",
        ],
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    output = result.stdout + result.stderr
    for domain in ("prok", "euk"):
        coverage_path = f"results/{sample}/coverm/{domain}_genome_coverage.tsv"
        assert (coverage_path in output) is (domain in expected_domains)
    assert f"results/{sample}/coverm/unknown_genome_coverage.tsv" not in output


def test_mandatory_coverm_validates_configured_input_mode_fields(tmp_path):
    config_path, _ = _coverm_config(tmp_path, ["prok"])
    config = json.loads(config_path.read_text(encoding="utf-8"))
    config["coverm"]["input_mode"] = "coupled"
    config_path.write_text(json.dumps(config), encoding="utf-8")

    result = subprocess.run(
        [
            "snakemake",
            "--snakefile",
            str(SNAKEFILE),
            "--configfile",
            str(config_path),
            "--directory",
            str(tmp_path),
            "--cores",
            "1",
            "--nolock",
            "--dry-run",
        ],
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
    assert "requires non-empty 'reads_1' for CoverM input_mode 'coupled'" in (
        result.stdout + result.stderr
    )


@pytest.mark.parametrize(
    ("entry_name", "entry_is_directory"),
    [
        (".snakemake_timestamp", False),
        ("bin.fna", False),
        ("nested.fa", True),
    ],
)
def test_coverm_skips_without_matching_regular_bin(
    tmp_path, monkeypatch, entry_name, entry_is_directory
):
    bins = tmp_path / "empty_bins"
    bins.mkdir()
    entry = bins / entry_name
    if entry_is_directory:
        entry.mkdir()
    else:
        entry.touch()
    output = tmp_path / "coverage.tsv"
    log = tmp_path / "logs" / "coverm.log"
    bam = tmp_path / "mapping.bam"
    bam.touch()
    fake_snakemake = SimpleNamespace(
        input=SimpleNamespace(bins=str(bins), bam=str(bam)),
        output=[str(output)],
        log=[str(log)],
        params=SimpleNamespace(
            input_mode="bam",
            mapper="minimap2-sr",
            methods=["mean"],
            extension=".fa",
            extra="",
        ),
        threads=1,
        log_fmt_shell=lambda **kwargs: f"> {log} 2>&1",
    )
    commands = []

    shell_module = importlib.import_module("snakemake.shell")
    monkeypatch.setattr(shell_module, "shell", commands.append)
    runpy.run_path(
        str(SCRIPTS / "coverm.py"),
        run_name="__main__",
        init_globals={"snakemake": fake_snakemake},
    )

    assert commands == []
    assert output.read_bytes() == b""
    assert log.read_text(encoding="utf-8") == (
        "CoverM skipped: no direct-child regular files matching genome FASTA "
        f"extension '.fa' in {bins}.\n"
    )


@pytest.mark.parametrize("native_output", [None, b"", b"Genome\n"])
def test_coverm_matching_bin_preserves_native_output_handling(
    tmp_path, monkeypatch, native_output
):
    bins = tmp_path / "bins"
    bins.mkdir()
    (bins / "bin.fa").write_text(">contig\nACGT\n", encoding="utf-8")
    output = tmp_path / "coverage.tsv"
    log = tmp_path / "coverm.log"
    bam = tmp_path / "mapping.bam"
    bam.touch()
    fake_snakemake = SimpleNamespace(
        input=SimpleNamespace(bins=str(bins), bam=str(bam)),
        output=[str(output)],
        log=[str(log)],
        params=SimpleNamespace(
            input_mode="bam",
            mapper="minimap2-sr",
            methods=["mean"],
            extension=".fa",
            extra="",
        ),
        threads=1,
        log_fmt_shell=lambda **kwargs: "> coverm.log 2>&1",
    )
    commands = []

    def fake_shell(command):
        commands.append(command)
        if native_output is not None:
            output.write_bytes(native_output)

    shell_module = importlib.import_module("snakemake.shell")
    monkeypatch.setattr(shell_module, "shell", fake_shell)
    runpy.run_path(
        str(SCRIPTS / "coverm.py"),
        run_name="__main__",
        init_globals={"snakemake": fake_snakemake},
    )

    assert len(commands) == 1
    assert commands[0].startswith("coverm genome ")
    if native_output is None:
        assert not output.exists()
    else:
        assert output.read_bytes() == native_output
    assert not log.exists()


def test_job_scripts_share_only_binning_specific_contracts():
    assert (
        prepare_refinement_inputs.ASSIGNMENT_COLUMNS
        is meew_binning_common.ASSIGNMENT_COLUMNS
    )
    assert run_das_tool.ASSIGNMENT_COLUMNS is meew_binning_common.ASSIGNMENT_COLUMNS
    assert (
        finalize_refined_domains.ASSIGNMENT_COLUMNS
        is meew_binning_common.ASSIGNMENT_COLUMNS
    )

    executable_modules = {
        "coverm",
        "finalize_refined_domains",
        "prepare_refinement_inputs",
        "run_das_tool",
        "semibin2",
        "separate_bins_by_domain",
    }
    for module in executable_modules:
        tree = ast.parse((SCRIPTS / f"{module}.py").read_text(encoding="utf-8"))
        imported_modules = {
            node.module.split(".", 1)[0]
            for node in ast.walk(tree)
            if isinstance(node, ast.ImportFrom) and node.module
        }
        imported_modules.update(
            alias.name.split(".", 1)[0]
            for node in ast.walk(tree)
            if isinstance(node, ast.Import)
            for alias in node.names
        )
        assert not imported_modules.intersection(executable_modules)
        assert "common" not in imported_modules

    assert not (SCRIPTS / "coverm_empty_outputs.py").exists()
    assert not (SCRIPTS / "meew_binning_identifiers.py").exists()
    assert not (SCRIPTS / "meew_binning_paths.py").exists()
