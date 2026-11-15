import csv
import os
import subprocess
import textwrap
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
TIARA_ENV = ROOT / "workflow" / "envs" / "tiara.yaml"
DOMAIN_RULES = ROOT / "workflow" / "rules" / "domain_separation.smk"

EXPECTED_DEPENDENCIES = [
    "python=3.8",
    "tiara=1.0.3=pyhd8ed1ab_0",
    "numpy=1.19.5=py38h8246c76_3",
    "pytorch=1.7.1=cpu_py38h36eccb8_2",
    "mkl=2020.4=h726a3e6_304",
    "libblas=3.9.0=8_mkl",
    "libcblas=3.9.0=8_mkl",
    "liblapack=3.9.0=8_mkl",
]


def test_tiara_environment_pins_compatible_cpu_mkl_stack():
    assert yaml.safe_load(TIARA_ENV.read_text(encoding="utf-8")) == {
        "channels": ["conda-forge", "bioconda", "nodefaults"],
        "dependencies": EXPECTED_DEPENDENCIES,
    }

    builds = {
        dependency.split("=", 1)[0]: dependency.rsplit("=", 1)[-1]
        for dependency in EXPECTED_DEPENDENCIES
        if dependency.startswith(("libblas=", "libcblas=", "liblapack="))
    }
    assert builds == {
        "libblas": "8_mkl",
        "libcblas": "8_mkl",
        "liblapack": "8_mkl",
    }


def test_tiara_rule_keeps_environment_threads_and_native_cli():
    tiara_rule = DOMAIN_RULES.read_text(encoding="utf-8").split(
        "\n\nif REFINEMENT_ENABLED:", 1
    )[0]
    assert 'conda:\n        "../envs/tiara.yaml"' in tiara_rule
    assert "threads: 4" in tiara_rule
    assert (
        "tiara -i {input.contigs:q} -o {output.tsv:q} --threads {threads} "
        "{params.extra} >{log:q} 2>&1"
    ) in tiara_rule


def _supplied_environment_prefix():
    value = os.environ.get("MEEW_TIARA_ENV_PREFIX")
    if not value:
        pytest.skip("set MEEW_TIARA_ENV_PREFIX to test a solved Tiara environment")
    prefix = Path(value)
    for executable in (prefix / "bin" / "python", prefix / "bin" / "tiara"):
        if not executable.is_file():
            pytest.fail(f"MEEW_TIARA_ENV_PREFIX is missing executable: {executable}")
    return prefix


def _smoke_environment(tmp_path):
    env = os.environ.copy()
    env.update(
        {
            "PYTHONDONTWRITEBYTECODE": "1",
            "TMPDIR": str(tmp_path),
            "XDG_CACHE_HOME": str(tmp_path / "cache"),
            "NUMBA_CACHE_DIR": str(tmp_path / "numba"),
        }
    )
    return env


def test_supplied_tiara_environment_native_imports_and_linear_algebra(tmp_path):
    prefix = _supplied_environment_prefix()
    result = subprocess.run(
        [
            str(prefix / "bin" / "python"),
            "-c",
            textwrap.dedent("""
                import numpy as np
                import scipy
                from scipy.linalg.blas import dgemm
                import torch
                import tiara

                assert np.__version__ == "1.19.5"
                assert torch.__version__ == "1.7.1"
                matrix = np.array([[4.0, 1.0], [2.0, 3.0]])
                vector = np.array([9.0, 13.0])
                assert np.allclose(np.linalg.solve(matrix, vector), [1.4, 3.4])
                assert np.allclose(dgemm(1.0, matrix, np.eye(2)), matrix)
                tensor = torch.tensor(matrix)
                assert torch.allclose(torch.mm(tensor, torch.eye(2)), tensor)
                print(np.__version__, scipy.__version__, torch.__version__)
                """),
        ],
        check=False,
        capture_output=True,
        text=True,
        env=_smoke_environment(tmp_path),
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_supplied_tiara_environment_native_prediction_matches_consumer_schema(tmp_path):
    prefix = _supplied_environment_prefix()
    fasta = tmp_path / "representative.fasta"
    output = tmp_path / "tiara.tsv"
    records = {
        "balanced_sequence": ("ACGTGCAATGCC" * 500),
        "gc_rich_sequence": ("GCGCGGCCGCGT" * 500),
        "at_rich_sequence": ("ATATTAATATGC" * 500),
    }
    fasta.write_text(
        "".join(f">{name}\n{sequence}\n" for name, sequence in records.items()),
        encoding="utf-8",
    )

    result = subprocess.run(
        [
            str(prefix / "bin" / "tiara"),
            "-i",
            str(fasta),
            "-o",
            str(output),
            "--threads",
            "4",
        ],
        check=False,
        capture_output=True,
        text=True,
        env=_smoke_environment(tmp_path),
    )
    assert result.returncode == 0, result.stdout + result.stderr

    with output.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    assert rows[0] == ["sequence_id", "class_fst_stage", "class_snd_stage"]
    assert [row[0] for row in rows[1:]] == list(records)
    assert all(len(row) == 3 for row in rows[1:])
    assert all(row[1] for row in rows[1:])
    assert all(row[2] for row in rows[1:])

    scripts = ROOT / "workflow" / "scripts"
    consumer = subprocess.run(
        [
            str(prefix / "bin" / "python"),
            "-c",
            textwrap.dedent("""
                import sys
                from pathlib import Path

                sys.path.insert(0, sys.argv[1])
                from meew_binning_common import read_tiara_labels

                labels = read_tiara_labels(Path(sys.argv[2]), "unknown")
                assert set(labels) == {
                    "balanced_sequence",
                    "gc_rich_sequence",
                    "at_rich_sequence",
                }
                assert set(labels.values()) <= {"prok", "euk", "unknown"}
                print(labels)
                """),
            str(scripts),
            str(output),
        ],
        check=False,
        capture_output=True,
        text=True,
        env=_smoke_environment(tmp_path),
    )
    assert consumer.returncode == 0, consumer.stdout + consumer.stderr
