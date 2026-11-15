import copy
import csv
import runpy
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml
from jsonschema import Draft202012Validator

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from meew_binning_common import (  # noqa: E402
    ASSIGNMENT_COLUMNS,
    generated_mag_identifier,
    normalize_domain,
    select_domain,
    validate_generated_mag_identifier,
)
from run_das_tool import assignment_counts  # noqa: E402


def _write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _read_tsv(path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _domain_outputs(root):
    return SimpleNamespace(
        prok_bins=str(root / "prok_bins"),
        euk_bins=str(root / "euk_bins"),
        unknown_bins=str(root / "unknown_bins"),
        summary=str(root / "domain_assignment.tsv"),
        annotation_samples=str(root / "annotation_samples.tsv"),
        rename_map=str(root / "bin_rename_map.tsv"),
    )


def _run_script(name, fake_snakemake):
    return runpy.run_path(
        str(SCRIPTS / name),
        run_name="__main__",
        init_globals={"snakemake": fake_snakemake},
    )


def test_refinement_script_ignores_conflicting_common_package(tmp_path):
    conflict_root = tmp_path / "conflict"
    _write(
        conflict_root / "common" / "__init__.py",
        "ASSIGNMENT_COLUMNS = ['conflicting']\n"
        "def wrong_helper(*args, **kwargs):\n"
        "    raise AssertionError('conflicting common helper was used')\n"
        "direct_child_fastas = wrong_helper\n"
        "read_tiara_labels = wrong_helper\n"
        "select_domain = wrong_helper\n"
        "serialize_provenance_path = wrong_helper\n"
        "unique_fasta_lengths = wrong_helper\n",
    )
    cached_script = tmp_path / "cache" / "prepare_refinement_inputs.py"
    _write(
        cached_script,
        (SCRIPTS / "prepare_refinement_inputs.py").read_text(encoding="utf-8"),
    )
    fasta = tmp_path / "one.fa"
    _write(fasta, ">contig\nACGT\n")

    result = subprocess.run(
        [
            sys.executable,
            "-c",
            "\n".join(
                [
                    "import importlib.util",
                    "import sys",
                    "from pathlib import Path",
                    "from types import SimpleNamespace",
                    f"sys.path.insert(0, {str(conflict_root)!r})",
                    "import common",
                    "assert common.ASSIGNMENT_COLUMNS == ['conflicting']",
                    f"script = Path({str(cached_script)!r})",
                    "spec = importlib.util.spec_from_file_location("
                    "'cached_prepare_refinement_inputs', script)",
                    "module = importlib.util.module_from_spec(spec)",
                    f"module.snakemake = SimpleNamespace(scriptdir={str(SCRIPTS)!r})",
                    "spec.loader.exec_module(module)",
                    "import meew_binning_common as binning_common",
                    "assert module.ASSIGNMENT_COLUMNS is "
                    "binning_common.ASSIGNMENT_COLUMNS",
                    "assert module.direct_child_fastas is "
                    "binning_common.direct_child_fastas",
                    "assert module.unique_fasta_lengths is "
                    "binning_common.unique_fasta_lengths",
                    f"fasta = Path({str(fasta)!r})",
                    "assert module.unique_fasta_lengths(fasta) == {'contig': 4}",
                    "assert module.direct_child_fastas(fasta.parent) == [fasta]",
                ]
            ),
        ],
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stdout + result.stderr


def test_schema_and_shipped_configs_default_to_unknown():
    schema = yaml.safe_load(
        (ROOT / "workflow" / "schemas" / "config.schema.yaml").read_text(
            encoding="utf-8"
        )
    )
    fallback_schema = schema["properties"]["domain_separation"]["properties"][
        "fallback_domain"
    ]
    assert fallback_schema["enum"] == ["prok", "euk", "unknown"]
    assert fallback_schema["default"] == "unknown"

    validator = Draft202012Validator(schema)
    for relative_path in (
        "config/config.yaml",
        ".test/config/config.yaml",
        ".test/config/config.single.yaml",
    ):
        config = yaml.safe_load((ROOT / relative_path).read_text(encoding="utf-8"))
        assert config["domain_separation"]["fallback_domain"] == "unknown"
        validator.validate(config)

    config = yaml.safe_load((ROOT / "config" / "config.yaml").read_text())
    for compatible_fallback in ("prok", "euk"):
        compatible_config = copy.deepcopy(config)
        compatible_config["domain_separation"]["fallback_domain"] = compatible_fallback
        validator.validate(compatible_config)


@pytest.mark.parametrize(
    ("bp_by_domain", "expected"),
    [
        ({}, "unknown"),
        ({"prok": 10, "euk": 10}, "unknown"),
        ({"prok": 1, "euk": 0, "unknown": 100}, "prok"),
        ({"prok": 0, "euk": 1, "unknown": 100}, "euk"),
    ],
)
def test_unknown_fallback_keeps_ties_safe_and_positive_calls(bp_by_domain, expected):
    assert select_domain(bp_by_domain, "unknown") == expected


def test_tiara_unknown_labels_and_explicit_legacy_fallbacks_are_supported():
    assert normalize_domain("unknown", "", "unknown") == "unknown"
    assert normalize_domain("unknown", "", "prok") == "prok"
    assert normalize_domain("unknown", "", "euk") == "euk"
    assert select_domain({}, "prok") == "prok"
    assert select_domain({"prok": 4, "euk": 4}, "euk") == "euk"


def test_generated_unknown_identifier_preserves_existing_validation():
    identifier = generated_mag_identifier("sample._-9", "unknown", 1)
    assert identifier == "sample._-9_unknown_bin_001"
    assert validate_generated_mag_identifier(identifier) == identifier

    with pytest.raises(ValueError, match="Invalid user identifier"):
        generated_mag_identifier("bad/sample", "unknown", 1)
    with pytest.raises(ValueError, match="Invalid generated MAG number"):
        generated_mag_identifier("sample", "unknown", 0)


def test_direct_routing_materializes_unknown_with_provenance_but_no_annotation(
    tmp_path, monkeypatch
):
    bins_dir = tmp_path / "source_bins"
    _write(bins_dir / "a_prok.fa", ">prok_contig\nAAAA\n")
    _write(bins_dir / "b_euk.fa", ">euk_contig\nAAAA\n")
    _write(
        bins_dir / "c_tied.fa",
        ">tie_prok\nAAAA\n>tie_euk\nTTTT\n",
    )
    _write(bins_dir / "d_unknown.fa", ">unknown_contig\nAAAA\n")
    _write(bins_dir / "e_no_evidence.fa", ">unlabelled_contig\nAAAA\n")
    tiara = tmp_path / "tiara.tsv"
    _write(
        tiara,
        "sequence_id\tclass_fst_stage\tclass_snd_stage\n"
        "prok_contig\tprokarya\t\n"
        "euk_contig\teukarya\t\n"
        "tie_prok\tprokarya\t\n"
        "tie_euk\teukarya\t\n"
        "unknown_contig\tunknown\t\n",
    )
    output_root = tmp_path / "domain"
    stale_unknown = output_root / "unknown_bins" / "stale.fa"
    _write(stale_unknown, ">stale\nA\n")
    fake_snakemake = SimpleNamespace(
        input=SimpleNamespace(tiara_tsv=str(tiara)),
        output=_domain_outputs(output_root),
        params=SimpleNamespace(bins_dir=str(bins_dir), fallback_domain="unknown"),
        wildcards=SimpleNamespace(sample="sample"),
    )
    monkeypatch.chdir(tmp_path)

    _run_script("separate_bins_by_domain.py", fake_snakemake)

    assert sorted(path.name for path in (output_root / "prok_bins").glob("*.fa")) == [
        "sample_prok_bin_001.fa"
    ]
    assert sorted(path.name for path in (output_root / "euk_bins").glob("*.fa")) == [
        "sample_euk_bin_001.fa"
    ]
    assert sorted(
        path.name for path in (output_root / "unknown_bins").glob("*.fa")
    ) == [
        "sample_unknown_bin_001.fa",
        "sample_unknown_bin_002.fa",
        "sample_unknown_bin_003.fa",
    ]
    assert not stale_unknown.exists()

    summary_rows = _read_tsv(output_root / "domain_assignment.tsv")
    assert [row["assigned_domain"] for row in summary_rows] == [
        "prok",
        "euk",
        "unknown",
        "unknown",
        "unknown",
    ]
    assert [row["standardized_bin"] for row in summary_rows[2:]] == [
        "sample_unknown_bin_001.fa",
        "sample_unknown_bin_002.fa",
        "sample_unknown_bin_003.fa",
    ]
    rename_rows = _read_tsv(output_root / "bin_rename_map.tsv")
    assert [row["domain"] for row in rename_rows] == [
        "prok",
        "euk",
        "unknown",
        "unknown",
        "unknown",
    ]
    annotation_rows = _read_tsv(output_root / "annotation_samples.tsv")
    assert {row["domain"] for row in annotation_rows} == {"prok", "euk"}
    assert all("unknown" not in row["sample"] for row in annotation_rows)


def _preparation_snakemake(tmp_path, bins_dirs, tools):
    output_root = tmp_path / "refinement"
    return SimpleNamespace(
        input=SimpleNamespace(
            bins_dirs=[str(path) for path in bins_dirs],
            assembly=str(tmp_path / "assembly.fa"),
            tiara_tsv=str(tmp_path / "tiara.tsv"),
        ),
        output=SimpleNamespace(
            mappings=[
                str(output_root / "das_tool" / "inputs" / f"{tool}_prok_bins.tsv")
                for tool in tools
            ],
            assignments=str(output_root / "bin_domain_assignments.tsv"),
        ),
        log=[str(tmp_path / "prepare.log")],
        params=SimpleNamespace(tools=tools, fallback_domain="unknown"),
        wildcards=SimpleNamespace(sample="sample"),
    )


def test_refinement_preserves_source_unknown_and_excludes_all_unknown_from_dastool(
    tmp_path, monkeypatch
):
    comebin_bins = tmp_path / "binning" / "comebin" / "bins"
    metabat_bins = tmp_path / "binning" / "metabat2" / "bins"
    _write(comebin_bins / "a_prok.fa", ">prok_contig\nAAAA\n")
    euk_bytes = b">euk_contig\r\nAAAA\r\n"
    unknown_bytes = b">unknown_contig\r\nAAAA\r\n"
    (comebin_bins / "b_euk.fa").parent.mkdir(parents=True, exist_ok=True)
    (comebin_bins / "b_euk.fa").write_bytes(euk_bytes)
    (comebin_bins / "c_unknown.fa").write_bytes(unknown_bytes)
    _write(metabat_bins / "a_other_unknown.fa", ">other_unknown\nAAAA\n")
    _write(metabat_bins / "b_other_prok.fa", ">other_prok\nAAAA\n")
    _write(
        tmp_path / "assembly.fa",
        ">prok_contig\nAAAA\n"
        ">euk_contig\nAAAA\n"
        ">unknown_contig\nAAAA\n"
        ">other_unknown\nAAAA\n"
        ">other_prok\nAAAA\n",
    )
    _write(
        tmp_path / "tiara.tsv",
        "sequence_id\tclass_fst_stage\tclass_snd_stage\n"
        "prok_contig\tprokarya\t\n"
        "euk_contig\teukarya\t\n"
        "unknown_contig\tunknown\t\n"
        "other_unknown\tunknown\t\n"
        "other_prok\tprokarya\t\n",
    )
    tools = ["comebin", "metabat2"]
    preparation = _preparation_snakemake(tmp_path, [comebin_bins, metabat_bins], tools)
    monkeypatch.chdir(tmp_path)

    _run_script("prepare_refinement_inputs.py", preparation)

    comebin_mapping = Path(preparation.output.mappings[0]).read_text(encoding="utf-8")
    metabat_mapping = Path(preparation.output.mappings[1]).read_text(encoding="utf-8")
    assert comebin_mapping == "prok_contig\tsample__comebin__bin_000001\n"
    assert metabat_mapping == "other_prok\tsample__metabat2__bin_000002\n"
    assert "unknown" not in comebin_mapping
    assert "unknown" not in metabat_mapping

    assignment_rows = _read_tsv(Path(preparation.output.assignments))
    assert [row["assigned_domain"] for row in assignment_rows] == [
        "prok",
        "euk",
        "unknown",
        "unknown",
        "prok",
    ]
    counts, expected_prokaryotic_ids = assignment_counts(
        Path(preparation.output.assignments), tools
    )
    assert counts["comebin"] == {
        "source_bins": 3,
        "prokaryotic_bins": 1,
        "eukaryotic_bins": 1,
        "unknown_bins": 1,
    }
    assert counts["metabat2"]["unknown_bins"] == 1
    assert expected_prokaryotic_ids == {
        "comebin": {"sample__comebin__bin_000001"},
        "metabat2": {"sample__metabat2__bin_000002"},
    }
    preparation_log = Path(preparation.log[0]).read_text(encoding="utf-8")
    assert "comebin.unknown_bins=1" in preparation_log
    assert "metabat2.unknown_bins=1" in preparation_log

    refined_prok = tmp_path / "dastool_bins"
    selected_bytes = b">prok_contig\r\nAAAA\r\n"
    refined_prok.mkdir()
    (refined_prok / "selected.fa").write_bytes(selected_bytes)
    output_root = tmp_path / "domain"
    for domain in ("prok", "euk", "unknown"):
        _write(output_root / f"{domain}_bins" / "stale.fa", ">stale\nA\n")
    finalization = SimpleNamespace(
        input=SimpleNamespace(
            prok_bins=str(refined_prok),
            source_bins=[str(comebin_bins), str(metabat_bins)],
            bin_domain_assignments=preparation.output.assignments,
            tiara_tsv=str(tmp_path / "tiara.tsv"),
        ),
        output=_domain_outputs(output_root),
        log=[str(tmp_path / "finalize.log")],
        params=SimpleNamespace(
            tools=tools,
            euk_source_tool="comebin",
            fallback_domain="unknown",
        ),
        wildcards=SimpleNamespace(sample="sample"),
    )

    _run_script("finalize_refined_domains.py", finalization)

    assert (output_root / "prok_bins" / "sample_prok_bin_001.fa").read_bytes() == (
        selected_bytes
    )
    assert (output_root / "euk_bins" / "sample_euk_bin_001.fa").read_bytes() == (
        euk_bytes
    )
    assert (
        output_root / "unknown_bins" / "sample_unknown_bin_001.fa"
    ).read_bytes() == unknown_bytes
    assert not any(path.name == "stale.fa" for path in output_root.rglob("*.fa"))

    summary_rows = _read_tsv(output_root / "domain_assignment.tsv")
    assert [row["original_bin"] for row in summary_rows] == [
        "selected.fa",
        "b_euk.fa",
        "c_unknown.fa",
    ]
    assert [row["assigned_domain"] for row in summary_rows] == [
        "prok",
        "euk",
        "unknown",
    ]
    assert "a_other_unknown.fa" not in {row["original_bin"] for row in summary_rows}
    rename_rows = _read_tsv(output_root / "bin_rename_map.tsv")
    assert [row["domain"] for row in rename_rows] == ["prok", "euk", "unknown"]
    annotation_rows = _read_tsv(output_root / "annotation_samples.tsv")
    assert [row["domain"] for row in annotation_rows] == ["prok", "euk"]
    assert Path(finalization.log[0]).read_text(encoding="utf-8") == (
        "final_prok_bins=1\nfinal_euk_bins=1\nfinal_unknown_bins=1\n"
    )


def test_das_tool_adapter_validates_and_reports_unknown_without_scientific_input(
    tmp_path, monkeypatch
):
    assembly = tmp_path / "assembly.fa"
    _write(assembly, ">unknown_contig\nAAAA\n")
    mapping = tmp_path / "inputs" / "comebin_prok_bins.tsv"
    _write(mapping, "")
    assignments = tmp_path / "bin_domain_assignments.tsv"
    assignments.parent.mkdir(parents=True, exist_ok=True)
    with assignments.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=ASSIGNMENT_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerow(
            {
                "tool": "comebin",
                "source_bin": "unknown.fa",
                "source_path": "source/unknown.fa",
                "normalized_bin_id": "sample__comebin__bin_000001",
                "assigned_domain": "unknown",
                "prok_bp": 0,
                "euk_bp": 0,
                "labeled_bp": 4,
                "unlabeled_bp": 0,
                "total_bp": 4,
                "n_contigs": 1,
                "n_labeled_contigs": 1,
                "fallback_domain": "unknown",
            }
        )

    output_root = tmp_path / "das_tool"
    fake_snakemake = SimpleNamespace(
        input=SimpleNamespace(
            assembly=str(assembly),
            contig2bin=[str(mapping)],
            bin_domain_assignments=str(assignments),
        ),
        output=SimpleNamespace(
            bins=str(output_root / "bins"),
            summary=str(output_root / "DASTool_summary.tsv"),
            scores=str(output_root / "all_bins.eval.tsv"),
            contig2bin=str(output_root / "DASTool_contig2bin.tsv"),
            provenance=str(output_root / "input_provenance.tsv"),
            run_status=str(output_root / "run_status.tsv"),
            tool_log=str(output_root / "DAS_Tool.log"),
        ),
        log=[str(tmp_path / "run.log")],
        params=SimpleNamespace(
            labels=["comebin"],
            expected_version="1.1.7",
            search_engine="diamond",
            score_threshold=0.5,
        ),
        threads=2,
    )
    commands = []

    def fake_run(command, **kwargs):
        commands.append(command)
        assert kwargs == {"check": True, "capture_output": True, "text": True}
        return SimpleNamespace(stdout="DAS Tool 1.1.7\n", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)

    _run_script("run_das_tool.py", fake_snakemake)

    assert commands == [["DAS_Tool", "--version"]]
    provenance_rows = _read_tsv(Path(fake_snakemake.output.provenance))
    assert provenance_rows == [
        {
            "label": "comebin",
            "contig2bin": str(mapping),
            "bin_domain_assignments": str(assignments),
            "source_bins": "1",
            "prokaryotic_bins": "0",
            "eukaryotic_bins": "0",
            "unknown_bins": "1",
            "prokaryotic_contigs": "0",
            "mapped_prokaryotic_bins": "0",
        }
    ]
    assert "unknownaryotic_bins" not in Path(
        fake_snakemake.output.provenance
    ).read_text(encoding="utf-8")
    assert "no_prokaryotic_inputs" in Path(fake_snakemake.output.run_status).read_text(
        encoding="utf-8"
    )
