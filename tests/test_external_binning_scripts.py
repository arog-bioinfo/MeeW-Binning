import json
import os
import re
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
RULES = ROOT / "workflow" / "rules" / "binning.smk"
COMEBIN_ATOMIC = ROOT / "workflow" / "scripts" / "run_comebin_atomic.sh"
METABINNER_PROFILES = ROOT / "workflow" / "scripts" / "prepare_metabinner_profiles.sh"


def _executable(path, content):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(textwrap.dedent(content).lstrip(), encoding="utf-8")
    path.chmod(0o755)
    return path


def _comebin_case(tmp_path, *, mode="success", runner_exit=""):
    case_root = tmp_path / "COMEBin case with spaces"
    bin_dir = case_root / "mock bin"
    bin_dir.mkdir(parents=True)
    _executable(
        bin_dir / "ln",
        """
        #!/usr/bin/env bash
        exit 1
        """,
    )
    _executable(
        bin_dir / "run_comebin.sh",
        f"""
        #!{sys.executable}
        import json
        import os
        import sys
        from pathlib import Path

        arguments = sys.argv[1:]
        Path(os.environ["NATIVE_CAPTURE"]).write_text(
            json.dumps(arguments), encoding="utf-8"
        )
        assembly = Path(arguments[arguments.index("-a") + 1])
        assert assembly.read_text(encoding="utf-8") == ">contig\\nACGT\\n"
        assert not assembly.is_symlink()
        output = Path(arguments[arguments.index("-o") + 1])
        result_root = output / "comebin_res"
        result_root.mkdir(parents=True, exist_ok=True)
        mode = os.environ["COMEBIN_TEST_MODE"]
        if mode != "missing_result":
            (result_root / "comebin_res.tsv").write_text(
                "contig\\tbin\\n", encoding="utf-8"
            )
        if mode != "missing_bins":
            bins = result_root / "comebin_res_bins"
            bins.mkdir()
            (bins / "bin 1.fa").write_text(">contig\\nACGT\\n", encoding="utf-8")
        """,
    )
    checkm_runner = _executable(
        case_root / "source cache" / "fake checkm runner.py",
        """
        import json
        import os
        import subprocess
        import sys
        from pathlib import Path

        arguments = sys.argv[1:]
        Path(os.environ["CHECKM_CAPTURE"]).write_text(
            json.dumps(arguments), encoding="utf-8"
        )
        if os.environ.get("CHECKM_RUNNER_EXIT"):
            raise SystemExit(int(os.environ["CHECKM_RUNNER_EXIT"]))
        separator = arguments.index("--")
        raise SystemExit(
            subprocess.run(arguments[separator + 1 :], check=False).returncode
        )
        """,
    )
    helper = case_root / "source cache" / "checkm helper.py"
    helper.write_text("# fixture\n", encoding="utf-8")
    source_assembly = case_root / "input data" / "assembly file.fa"
    source_assembly.parent.mkdir(parents=True)
    source_assembly.write_text(">contig\nACGT\n", encoding="utf-8")
    bam = case_root / "input data" / "alignment file.bam"
    bam.write_text("fixture\n", encoding="utf-8")
    final_output = case_root / "published results" / "comebin output"
    final_output.mkdir(parents=True)
    (final_output / "old sentinel").write_text("old\n", encoding="utf-8")
    checkm_capture = case_root / "checkm arguments.json"
    native_capture = case_root / "native arguments.json"
    database = case_root / "database with spaces"

    environment = os.environ.copy()
    environment.update(
        {
            "PATH": str(bin_dir) + os.pathsep + environment["PATH"],
            "CHECKM_CAPTURE": str(checkm_capture),
            "NATIVE_CAPTURE": str(native_capture),
            "COMEBIN_TEST_MODE": mode,
            "CHECKM_RUNNER_EXIT": runner_exit,
        }
    )
    result = subprocess.run(
        [
            "bash",
            str(COMEBIN_ATOMIC),
            "--checkm-runner",
            str(checkm_runner),
            "--checkm-helper",
            str(helper),
            "--database",
            str(database),
            "--download-database",
            "True",
            "--source-assembly",
            str(source_assembly),
            "--bam",
            str(bam),
            "--bam-dir",
            str(bam.parent),
            "--final-output",
            str(final_output),
            "--threads",
            "3",
            "--",
            "--extra-flag",
            "value with spaces",
        ],
        check=False,
        capture_output=True,
        text=True,
        env=environment,
    )
    return {
        "result": result,
        "root": case_root,
        "source": source_assembly,
        "bam": bam,
        "bam_dir": bam.parent,
        "helper": helper,
        "database": database,
        "final": final_output,
        "checkm_capture": checkm_capture,
        "native_capture": native_capture,
    }


def _temporary_workspaces(case):
    final_output = case["final"]
    return list(final_output.parent.glob(f".{final_output.name}.tmp.*"))


def test_comebin_atomic_script_preserves_cached_arguments_and_publishes(tmp_path):
    case = _comebin_case(tmp_path)
    result = case["result"]
    assert result.returncode == 0, result.stdout + result.stderr
    assert "hard link unavailable; copying" in result.stdout
    assert not _temporary_workspaces(case)
    assert not (case["final"] / "old sentinel").exists()
    assert (case["final"] / "comebin_res" / "comebin_res.tsv").read_text(
        encoding="utf-8"
    ) == ("contig\tbin\n")
    assert (case["final"] / "comebin_res" / "comebin_res_bins" / "bin 1.fa").is_file()

    checkm_arguments = json.loads(case["checkm_capture"].read_text(encoding="utf-8"))
    protected_paths = [
        checkm_arguments[index + 1]
        for index, argument in enumerate(checkm_arguments)
        if argument == "--protected-path"
    ]
    assert protected_paths[:3] == [
        str(case["source"]),
        str(case["bam"]),
        str(case["final"]),
    ]
    workspace = Path(protected_paths[3])
    assert workspace.parent == case["final"].parent
    assert checkm_arguments[:6] == [
        "--helper-source",
        str(case["helper"]),
        "--database",
        str(case["database"]),
        "--download-database",
        "True",
    ]
    separator = checkm_arguments.index("--")
    native_arguments = checkm_arguments[separator + 1 :]
    assert native_arguments == [
        "run_comebin.sh",
        "-a",
        str(workspace / "assembly" / case["source"].name),
        "-p",
        str(case["bam_dir"]),
        "-t",
        "3",
        "-o",
        str(workspace / "output"),
        "--extra-flag",
        "value with spaces",
    ]
    assert json.loads(case["native_capture"].read_text(encoding="utf-8")) == (
        native_arguments[1:]
    )


def test_comebin_atomic_script_cleans_workspace_after_runner_failure(tmp_path):
    case = _comebin_case(tmp_path, runner_exit="23")
    assert case["result"].returncode == 23
    assert (case["final"] / "old sentinel").read_text(encoding="utf-8") == "old\n"
    assert not _temporary_workspaces(case)
    assert not case["native_capture"].exists()


@pytest.mark.parametrize(
    ("mode", "message"),
    [
        ("missing_result", "COMEBin did not produce comebin_res/comebin_res.tsv"),
        ("missing_bins", "COMEBin did not produce comebin_res/comebin_res_bins"),
    ],
)
def test_comebin_atomic_script_validates_outputs_and_cleans(tmp_path, mode, message):
    case = _comebin_case(tmp_path, mode=mode)
    assert case["result"].returncode == 1
    assert message in case["result"].stdout
    assert (case["final"] / "old sentinel").read_text(encoding="utf-8") == "old\n"
    assert not _temporary_workspaces(case)


def test_metabinner_profile_script_preserves_filtering_and_outputs(tmp_path):
    case_root = tmp_path / "MetaBinner case with spaces"
    bin_dir = case_root / "mock bin"
    bin_dir.mkdir(parents=True)
    filter_capture = case_root / "filter arguments.json"
    kmer_capture = case_root / "kmer arguments.json"
    _executable(
        bin_dir / "Filter_tooshort.py",
        f"""
        #!{sys.executable}
        import json
        import os
        import sys
        from pathlib import Path

        source = Path(sys.argv[1])
        cutoff = sys.argv[2]
        Path(os.environ["FILTER_CAPTURE"]).write_text(
            json.dumps(sys.argv[1:]), encoding="utf-8"
        )
        source.with_name(f"{{source.stem}}_{{cutoff}}{{source.suffix}}").write_text(
            ">keep\\nAAAA\\n", encoding="utf-8"
        )
        source.with_name(f"contig_length_filter{{cutoff}}.txt").write_text(
            "keep\\t4\\n", encoding="utf-8"
        )
        """,
    )
    _executable(
        bin_dir / "gen_kmer.py",
        f"""
        #!{sys.executable}
        import json
        import os
        import sys
        from pathlib import Path

        source = Path(sys.argv[1])
        Path(os.environ["KMER_CAPTURE"]).write_text(
            json.dumps(sys.argv[1:]), encoding="utf-8"
        )
        source.with_name("filtered_contigs_kmer_4_f0.csv").write_text(
            "contig,kmer\\nkeep,1\\n", encoding="utf-8"
        )
        """,
    )
    inputs = case_root / "input data"
    inputs.mkdir()
    contigs = inputs / "assembly file.fa"
    contigs.write_text(">keep\nAAAA\n>drop\nAA\n", encoding="utf-8")
    depth = inputs / "depth profile.tsv"
    depth.write_text(
        "contigName\tcontigLen\ttotalAvgDepth\tsample one\tsample two\n"
        "keep\t4\t7.5\t7\t8\n"
        "drop\t2\t2.5\t2\t3\n",
        encoding="utf-8",
    )
    output_dir = case_root / "output profiles"
    output_contigs = output_dir / "filtered_contigs.fa"
    output_coverage = output_dir / "coverage_profile.tsv"
    output_composition = output_dir / "tetranucleotide_profile.csv"
    environment = os.environ.copy()
    environment.update(
        {
            "PATH": str(bin_dir) + os.pathsep + environment["PATH"],
            "FILTER_CAPTURE": str(filter_capture),
            "KMER_CAPTURE": str(kmer_capture),
        }
    )

    result = subprocess.run(
        [
            "bash",
            str(METABINNER_PROFILES),
            "--contigs",
            str(contigs),
            "--depth",
            str(depth),
            "--output-dir",
            str(output_dir),
            "--output-contigs",
            str(output_contigs),
            "--output-coverage",
            str(output_coverage),
            "--output-composition",
            str(output_composition),
            "--min-contig-length",
            "1000",
        ],
        check=False,
        capture_output=True,
        text=True,
        env=environment,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert output_contigs.read_text(encoding="utf-8") == ">keep\nAAAA\n"
    assert output_coverage.read_text(encoding="utf-8") == (
        "contigName\tsample one\tsample two\nkeep\t7\t8\n"
    )
    assert output_composition.read_text(encoding="utf-8") == ("contig,kmer\nkeep,1\n")
    assert json.loads(filter_capture.read_text(encoding="utf-8")) == [
        str(output_dir / "metabinner_input.fa"),
        "1000",
    ]
    assert json.loads(kmer_capture.read_text(encoding="utf-8")) == [
        str(output_contigs),
        "0",
        "4",
    ]
    assert not (output_dir / "metabinner_input.fa").exists()
    assert not (output_dir / "contig_length_filter1000.txt").exists()


def test_rule_shell_blocks_are_brace_free_and_scripts_are_cached_inputs():
    source = RULES.read_text(encoding="utf-8")
    shell_blocks = re.findall(r'shell:\n\s+"""\n(.*?)\n\s+"""', source, re.DOTALL)
    assert shell_blocks
    assert all("{{" not in block and "}}" not in block for block in shell_blocks)
    assert 'atomic_script=workflow.source_path("../scripts/run_comebin_atomic.sh")' in (
        source
    )
    assert (
        'profile_script=workflow.source_path("../scripts/prepare_metabinner_profiles.sh")'
        in source
    )
    assert "command -v run_comebin.sh" not in "\n".join(shell_blocks)
    assert "${{" not in "\n".join(shell_blocks)
