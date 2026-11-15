from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
WORKFLOW = ROOT / ".github" / "workflows" / "main.yaml"


def _load_workflow():
    return yaml.load(WORKFLOW.read_text(encoding="utf-8"), Loader=yaml.BaseLoader)


def _named_step(job, name):
    return next(step for step in job["steps"] if step.get("name") == name)


def test_template_plus_ci_contract():
    workflow = _load_workflow()

    triggers = workflow["on"]
    assert set(triggers) == {"push", "pull_request", "workflow_dispatch"}
    assert triggers["push"]["branches"] == ["main"]
    assert triggers["pull_request"]["branches"] == ["main"]
    assert workflow["permissions"] == {"contents": "read"}
    assert workflow["concurrency"] == {
        "group": "${{ github.workflow }}-${{ github.ref }}",
        "cancel-in-progress": "true",
    }

    jobs = workflow["jobs"]
    assert set(jobs) == {"formatting", "linting", "python-tests", "testing"}
    assert "needs" not in jobs["formatting"]
    assert "needs" not in jobs["linting"]
    assert "needs" not in jobs["python-tests"]
    assert jobs["testing"]["needs"] == [
        "formatting",
        "linting",
        "python-tests",
    ]
    assert {name: job["timeout-minutes"] for name, job in jobs.items()} == {
        "formatting": "15",
        "linting": "15",
        "python-tests": "20",
        "testing": "20",
    }
    for job in jobs.values():
        checkout_refs = [
            step["uses"]
            for step in job["steps"]
            if step.get("uses", "").startswith("actions/checkout@")
        ]
        assert checkout_refs == ["actions/checkout@v6"]

    formatting = jobs["formatting"]
    assert formatting["steps"][0] == {
        "uses": "actions/checkout@v6",
        "with": {"fetch-depth": "0"},
    }
    formatter = _named_step(formatting, "Check formatting")
    assert formatter["uses"] == "super-linter/super-linter@v8"
    assert {
        key: value
        for key, value in formatter["env"].items()
        if key.startswith("VALIDATE_")
    } == {
        "VALIDATE_ALL_CODEBASE": "true",
        "VALIDATE_PYTHON_BLACK": "true",
        "VALIDATE_SHELL_SHFMT": "true",
        "VALIDATE_SNAKEMAKE_SNAKEFMT": "true",
        "VALIDATE_YAML_PRETTIER": "true",
    }

    lint = _named_step(jobs["linting"], "Lint workflow")
    assert lint["uses"] == "snakemake/snakemake-github-action@v2"
    assert lint["with"] == {
        "directory": ".",
        "snakefile": "workflow/Snakefile",
        "args": "--lint --configfile .test/config/config.single.yaml",
        "snakemake-version": "9.19.0",
    }

    python_job = jobs["python-tests"]
    setup_python = next(
        step
        for step in python_job["steps"]
        if step.get("uses") == "actions/setup-python@v6"
    )
    assert setup_python["with"]["python-version"] == "3.13"
    assert setup_python["with"]["cache"] == "pip"
    assert setup_python["with"]["cache-dependency-path"].splitlines() == [
        "tests/requirements.txt",
        "tests/requirements-workflow.txt",
    ]
    assert _named_step(python_job, "Install test dependencies")["run"] == (
        "python -m pip install -r tests/requirements-workflow.txt"
    )
    assert _named_step(python_job, "Run Python tests")["run"] == (
        "python -m pytest -q -ra tests"
    )

    suites = jobs["testing"]["strategy"]["matrix"]["suite"]
    assert suites == [
        {"name": "direct", "config": ".test/config/config.single.yaml"},
        {"name": "refinement", "config": ".test/config/config.refinement.yaml"},
    ]

    testing_steps = jobs["testing"]["steps"]
    assert testing_steps[0] == {"uses": "actions/checkout@v6"}
    assert len(testing_steps) == 2
    dry_run = _named_step(jobs["testing"], "Dry-run workflow")
    assert dry_run["uses"] == "snakemake/snakemake-github-action@v2"
    assert dry_run["with"]["directory"] == "."
    assert dry_run["with"]["snakefile"] == "workflow/Snakefile"
    assert dry_run["with"]["snakemake-version"] == "9.19.0"

    dag_args = " ".join(dry_run["with"]["args"].split())
    assert dag_args == (
        "--dry-run --printshellcmds --cores 2 "
        "--configfile ${{ matrix.suite.config }}"
    )
    for forbidden in (
        "--software-deployment-method",
        "--sdm",
        "conda",
        "create-env",
        "--report",
        "upload-artifact",
    ):
        assert forbidden not in dag_args.lower()


def test_ci_dependencies_are_pinned():
    requirements = (ROOT / "tests" / "requirements.txt").read_text(encoding="utf-8")
    workflow_requirements = (ROOT / "tests" / "requirements-workflow.txt").read_text(
        encoding="utf-8"
    )

    assert requirements.splitlines() == [
        "pytest==9.0.2",
        "PyYAML==6.0.3",
        "jsonschema==4.26.0",
        "pandas==2.3.3",
    ]
    assert workflow_requirements.splitlines() == [
        "-r requirements.txt",
        "snakemake==9.19.0",
    ]
