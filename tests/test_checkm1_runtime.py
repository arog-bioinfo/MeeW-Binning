import fcntl
import hashlib
import importlib.util
import json
import os
import signal
import subprocess
import sys
import threading
import time
from pathlib import Path

import pytest
import yaml
from jsonschema import Draft202012Validator

ROOT = Path(__file__).resolve().parents[1]
HELPER = ROOT / "workflow" / "scripts" / "common" / "checkm1.py"
RUNNER = ROOT / "workflow" / "scripts" / "run_comebin_with_checkm.py"
COMEBIN_ATOMIC = ROOT / "workflow" / "scripts" / "run_comebin_atomic.sh"
METABINNER_PROFILES = ROOT / "workflow" / "scripts" / "prepare_metabinner_profiles.sh"
RULES = ROOT / "workflow" / "rules" / "binning.smk"
SCHEMA = ROOT / "workflow" / "schemas" / "config.schema.yaml"


def _load_helper():
    spec = importlib.util.spec_from_file_location("meew_binning_checkm1_test", HELPER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


checkm1 = _load_helper()


def _sha256(payload):
    return hashlib.sha256(payload).hexdigest()


def _ready_database(root):
    payloads = {
        "hmms/checkm.hmm": (
            b"HMMER3/f [3.1b1 | May 2013]\n"
            b"NAME  fixture\n"
            b"ACC   PF00001.1\n"
            b"LENG  1\n"
            b"ALPH  amino\n"
            b"HMM          A        C\n"
            b"      1   0.0 0.0\n//\n"
        ),
        "hmms/checkm.hmm.ssi": (
            b"\xd3\xd3\xc9\xb3" + b"\x00" * 32 + b"checkm.hmm\x00" + b"index"
        ),
        "pfam/Pfam-A.hmm.dat": (
            b"# STOCKHOLM 1.0\n"
            b"#=GF ID   fixture\n"
            b"#=GF AC   PF00001.1\n"
            b"#=GF GA   1.00; 1.00;\n"
            b"#=GF TP   Family\n//\n"
        ),
    }
    manifest = ["##CheckM##\tData manifest created by ScreamingBackpack version 0.2.3"]
    for relative, payload in payloads.items():
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(payload)
        manifest.append("%s\t%s\t%d" % (relative, _sha256(payload), len(payload)))
    root.mkdir(parents=True, exist_ok=True)
    (root / ".dmanifest").write_text("\n".join(manifest) + "\n", encoding="utf-8")
    return root


def _write_marker(path, domain):
    markers = ["PF%05d.1" % value for value in range(1, 46)]
    expression = "[set([%s])]" % ", ".join(repr(marker) for marker in markers)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "# [Taxon Marker File]\n%s\t1\t0\t%s\t1\t%s\n" % (domain, domain, expression),
        encoding="utf-8",
    )


def _marker_root(root):
    _write_marker(root / "bacteria.ms", "Bacteria")
    _write_marker(root / "archaea.ms", "Archaea")
    return root


def _executable(path, content="#!/bin/sh\nexit 0\n"):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    path.chmod(0o755)
    return path


def _hook(path):
    return _executable(
        path,
        "#!/bin/bash\n"
        "# CHECKM_DATA_DIR CHECKM_DATA_PATH\n"
        "# %s\n# %s\n"
        "# sha256sum -c\n# checkm data setRoot\nexit 0\n"
        % (checkm1.OFFICIAL_ARCHIVE_URL, checkm1.OFFICIAL_ARCHIVE_SHA256),
    )


class _Forwarder:
    received = None

    def raise_if_cancelled(self):
        return None

    def run(self, command, environ=None):
        del command, environ
        return 0

    def run_capture(self, command, environ=None):
        del command, environ
        return 0, "HMMER3/f\nNAME fixture\n//\n", ""


def test_path_precedence_and_legacy_blank_defaults():
    environment = {
        "HOME": "/tmp/checkm home",
        "XDG_DATA_HOME": "/tmp/checkm xdg",
        "MEEW_RESOURCES": "/tmp/checkm ambient",
    }
    assert (
        checkm1.resolve_checkm1_database(
            "/tmp/exact db", "/tmp/module root", environment
        )
        == "/tmp/exact db"
    )
    assert (
        checkm1.resolve_checkm1_database("", "/tmp/module root", environment)
        == "/tmp/module root/checkm1_db"
    )
    assert (
        checkm1.resolve_checkm1_database("", "   ", environment)
        == "/tmp/checkm ambient/checkm1_db"
    )
    del environment["MEEW_RESOURCES"]
    assert (
        checkm1.resolve_checkm1_database("", "", environment)
        == "/tmp/checkm xdg/meew/checkm1_db"
    )
    del environment["XDG_DATA_HOME"]
    assert (
        checkm1.resolve_checkm1_database("", "", environment)
        == "/tmp/checkm home/resources/checkm1_db"
    )


def test_schema_accepts_omitted_legacy_keys_and_ships_auto_enabled_default():
    schema = yaml.safe_load(SCHEMA.read_text(encoding="utf-8"))
    validator = Draft202012Validator(schema)
    shipped = yaml.safe_load(
        (ROOT / "config" / "config.yaml").read_text(encoding="utf-8")
    )
    validator.validate(shipped)
    assert shipped["resources_root"] == ""
    assert shipped["checkm1"] == {
        "database_path": "",
        "download_database": True,
    }

    legacy = dict(shipped)
    legacy.pop("resources_root")
    legacy.pop("checkm1")
    validator.validate(legacy)
    assert checkm1._parse_boolean("true") is True
    assert checkm1._parse_boolean("false") is False


def test_database_readiness_rejects_arbitrary_nonempty_and_corrupt_payloads(tmp_path):
    database = tmp_path / "database"
    for relative in checkm1.REQUIRED_DATABASE_FILES:
        path = database / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("ready\n", encoding="utf-8")
    with pytest.raises(checkm1.CheckM1Error, match="manifest header"):
        checkm1.validate_database(database)

    _ready_database(database)
    checkm1.validate_database(database)

    model = database / "hmms" / "checkm.hmm"
    model.write_bytes(model.read_bytes()[:-3])
    with pytest.raises(checkm1.CheckM1Error, match="size is|truncated"):
        checkm1.validate_database(database)

    _ready_database(database)
    index = database / "hmms" / "checkm.hmm.ssi"
    payload = bytearray(index.read_bytes())
    payload[0] = 0
    index.write_bytes(payload)
    manifest = database / ".dmanifest"
    manifest.write_text(
        manifest.read_text(encoding="utf-8").replace(
            _sha256(bytes([0xD3]) + bytes(payload[1:])), _sha256(bytes(payload))
        ),
        encoding="utf-8",
    )
    with pytest.raises(checkm1.CheckM1Error, match="SSI index"):
        checkm1.validate_database(database)


def test_marker_readiness_rejects_nonempty_malformed_files(tmp_path):
    marker_root = _marker_root(tmp_path / "checkm_ms")
    assert checkm1._marker_file_problems(marker_root / "bacteria.ms", "Bacteria") == []
    (marker_root / "bacteria.ms").write_text("ready\n", encoding="utf-8")
    assert (
        "not a two-row"
        in checkm1._marker_file_problems(marker_root / "bacteria.ms", "Bacteria")[0]
    )


def test_nonempty_incomplete_or_corrupt_cache_is_never_replaced(tmp_path):
    destination = tmp_path / "database"
    destination.mkdir()
    sentinel = destination / "keep-me"
    sentinel.write_text("protected\n", encoding="utf-8")
    with pytest.raises(checkm1.CheckM1Error, match="nonempty incomplete/corrupt"):
        checkm1._destination_state(destination)
    assert sentinel.read_text(encoding="utf-8") == "protected\n"


def test_opt_out_rejects_missing_or_empty_destination(tmp_path):
    config_path = tmp_path / "checkm" / "DATA_CONFIG"
    config_path.parent.mkdir()
    snapshot = b'{"dataRoot": "/prior"}'
    config_path.write_bytes(snapshot)
    for destination in (tmp_path / "missing", tmp_path / "empty"):
        if destination.name == "empty":
            destination.mkdir()
        with pytest.raises(checkm1.CheckM1Error, match="download_database is false"):
            checkm1._prepare_database(
                destination,
                False,
                tmp_path / "hook",
                tmp_path / "environment",
                config_path,
                snapshot,
                _Forwarder(),
                {},
            )


def test_native_staging_publishes_atomically_and_keeps_registration_isolated(
    tmp_path, monkeypatch
):
    destination = tmp_path / "resources" / "checkm_database"
    config_path = tmp_path / "checkm" / "DATA_CONFIG"
    config_path.parent.mkdir()
    snapshot = b'{"dataRoot": "/prior"}'
    config_path.write_bytes(snapshot)
    calls = []

    def successful_hook(
        hook,
        staged_payload,
        environment_prefix,
        observed_config,
        observed_snapshot,
        forwarder,
    ):
        del hook, environment_prefix, forwarder
        assert observed_config == config_path
        assert observed_snapshot == snapshot
        calls.append(staged_payload)
        _ready_database(staged_payload)

    monkeypatch.setattr(checkm1, "_run_native_hook", successful_hook)
    monkeypatch.setattr(checkm1, "_native_database_probe", lambda *args: None)
    checkm1._prepare_database(
        destination,
        True,
        tmp_path / "hook",
        tmp_path / "environment",
        config_path,
        snapshot,
        _Forwarder(),
        {},
    )
    assert calls
    checkm1.validate_database(destination)
    assert config_path.read_bytes() == snapshot
    assert not any(
        path.name.startswith(".checkm_database.meew-stage-")
        for path in destination.parent.iterdir()
    )


@pytest.mark.parametrize("failure", ["network", "checksum"])
def test_download_or_checksum_failure_leaves_no_partial_cache(
    tmp_path, monkeypatch, failure
):
    destination = tmp_path / (failure + "_database")
    config_path = tmp_path / "checkm" / "DATA_CONFIG"
    config_path.parent.mkdir()
    snapshot = b'{"dataRoot": "/prior"}'
    config_path.write_bytes(snapshot)

    def failed_hook(*args, **kwargs):
        del args, kwargs
        raise checkm1.CheckM1Error(failure)

    monkeypatch.setattr(checkm1, "_run_native_hook", failed_hook)
    with pytest.raises(checkm1.CheckM1Error, match=failure):
        checkm1._prepare_database(
            destination,
            True,
            tmp_path / "hook",
            tmp_path / "environment",
            config_path,
            snapshot,
            _Forwarder(),
            {},
        )
    assert config_path.read_bytes() == snapshot
    assert not destination.exists()


def test_publication_failure_preserves_existing_state(tmp_path, monkeypatch):
    destination = tmp_path / "database"
    config_path = tmp_path / "checkm" / "DATA_CONFIG"
    config_path.parent.mkdir()
    snapshot = b'{"dataRoot": "/prior"}'
    config_path.write_bytes(snapshot)

    def successful_hook(hook, staged_payload, *args):
        del hook, args
        _ready_database(staged_payload)

    original_replace = checkm1.os.replace

    def denied(source, target):
        if Path(target) == destination:
            raise OSError("publication denied")
        return original_replace(source, target)

    monkeypatch.setattr(checkm1, "_run_native_hook", successful_hook)
    monkeypatch.setattr(checkm1, "_native_database_probe", lambda *args: None)
    monkeypatch.setattr(checkm1.os, "replace", denied)
    with pytest.raises(checkm1.CheckM1Error, match="atomically published"):
        checkm1._prepare_database(
            destination,
            True,
            tmp_path / "hook",
            tmp_path / "environment",
            config_path,
            snapshot,
            _Forwarder(),
            {},
        )
    assert not destination.exists()
    assert config_path.read_bytes() == snapshot


def test_valid_cache_reuse_skips_native_hook(tmp_path, monkeypatch):
    destination = _ready_database(tmp_path / "database")
    config_path = tmp_path / "checkm" / "DATA_CONFIG"
    config_path.parent.mkdir()
    snapshot = b'{"dataRoot": "/prior"}'
    config_path.write_bytes(snapshot)
    probes = []

    def forbidden(*args, **kwargs):
        raise AssertionError("native hook must not run for a valid cache")

    monkeypatch.setattr(checkm1, "_run_native_hook", forbidden)
    monkeypatch.setattr(
        checkm1, "_native_database_probe", lambda *args: probes.append(args[0])
    )
    checkm1._prepare_database(
        destination,
        True,
        tmp_path / "hook",
        tmp_path / "environment",
        config_path,
        snapshot,
        _Forwarder(),
        {},
    )
    assert probes == [destination]
    assert config_path.read_bytes() == snapshot


def test_native_hook_boundary_confines_prefix_archive_override_and_registration(
    tmp_path, monkeypatch
):
    staged = tmp_path / "payload with spaces"
    staged.mkdir()
    package = (
        tmp_path / "environment" / "lib" / "python3.7" / "site-packages" / "checkm"
    )
    package.mkdir(parents=True)
    config_path = package / "DATA_CONFIG"
    snapshot = b'{"dataRoot": "/prior"}'
    config_path.write_bytes(snapshot)
    monkeypatch.setenv("CHECKM_DATA_PATH", "/untrusted/local/archive")
    observed = {}

    class Forwarder(_Forwarder):
        def run(self, command, environ=None):
            observed["command"] = command
            observed["environment"] = dict(environ)
            isolated = Path(environ["PYTHONPATH"].split(os.pathsep)[0]) / "checkm"
            (isolated / "DATA_CONFIG").write_text(
                json.dumps({"dataRoot": environ["CHECKM_DATA_DIR"]}),
                encoding="utf-8",
            )
            _ready_database(Path(environ["CHECKM_DATA_DIR"]))
            return 0

    checkm1._run_native_hook(
        tmp_path / "installed_hook",
        staged,
        tmp_path / "environment",
        config_path,
        snapshot,
        Forwarder(),
    )
    assert observed["command"][0] == "bash"
    assert "CHECKM_DATA_PATH" not in observed["environment"]
    assert " " not in observed["environment"]["PREFIX"]
    assert " " not in observed["environment"]["CHECKM_DATA_DIR"]
    assert Path(observed["environment"]["HOME"]).name == "home"
    assert config_path.read_bytes() == snapshot
    checkm1.validate_database(staged)


def test_tool_specific_runner_and_marker_roots_do_not_treat_bash_as_install(tmp_path):
    prefix = tmp_path / "environment"
    bin_dir = prefix / "bin"
    _hook(bin_dir / ".checkm-genome-post-link.sh")
    for executable in set(
        checkm1.REQUIRED_SETUP_EXECUTABLES + checkm1.REQUIRED_SCORING_EXECUTABLES
    ):
        _executable(bin_dir / executable)
    environment = {"PATH": str(bin_dir)}

    comebin_runner = _executable(bin_dir / "run_comebin.sh")
    comebin_markers = _marker_root(bin_dir / "COMEBin" / "scripts" / "checkm_ms")
    comebin = checkm1.validate_native_installation(
        "comebin", [str(comebin_runner)], environment, prefix
    )
    assert comebin["marker_root"] == comebin_markers

    tool_root = tmp_path / "MetaBinner"
    metabinner_runner = _executable(tool_root / "run_metabinner.sh")
    (tool_root / "auxiliary").mkdir(parents=True)
    metabinner_markers = _marker_root(tool_root / "scripts" / "checkm_ms")
    metabinner = checkm1.validate_native_installation(
        "metabinner",
        [str(bin_dir / "bash"), str(metabinner_runner)],
        environment,
        prefix,
        str(tool_root),
    )
    assert metabinner["install_root"] == tool_root
    assert metabinner["marker_root"] == metabinner_markers
    assert metabinner["install_root"] != bin_dir

    (metabinner_markers / "archaea.ms").write_text("ready\n", encoding="utf-8")
    with pytest.raises(checkm1.CheckM1Error, match="marker readiness"):
        checkm1.validate_native_installation(
            "metabinner",
            [str(bin_dir / "bash"), str(metabinner_runner)],
            environment,
            prefix,
            str(tool_root),
        )


def test_missing_hmmfetch_and_hmmalign_are_reported(tmp_path):
    prefix = tmp_path / "environment"
    bin_dir = prefix / "bin"
    _hook(bin_dir / ".checkm-genome-post-link.sh")
    omitted = {"hmmfetch", "hmmalign"}
    for executable in (
        set(checkm1.REQUIRED_SETUP_EXECUTABLES + checkm1.REQUIRED_SCORING_EXECUTABLES)
        - omitted
    ):
        _executable(bin_dir / executable)
    runner = _executable(bin_dir / "run_comebin.sh")
    _marker_root(bin_dir / "COMEBin" / "scripts" / "checkm_ms")
    with pytest.raises(checkm1.CheckM1Error, match="hmmfetch, hmmalign"):
        checkm1.validate_native_installation(
            "comebin", [str(runner)], {"PATH": str(bin_dir)}, prefix
        )


def test_native_path_boundaries_reject_unsupported_whitespace_early():
    with pytest.raises(checkm1.CheckM1Error, match="Whitespace-containing"):
        checkm1.validate_native_paths("/tmp/checkm db", "comebin")
    with pytest.raises(checkm1.CheckM1Error, match="MetaBinner"):
        checkm1.validate_native_paths(
            "/tmp/checkm_db",
            "metabinner",
            "/tmp/meta binner",
            ["/tmp/input.fa"],
        )
    checkm1.validate_native_paths(
        "/tmp/checkm_db", "metabinner", "/tmp/metabinner", ["/tmp/input.fa"]
    )


def test_overlap_and_shared_environment_lock_identity(tmp_path):
    destination = tmp_path / "results" / "checkm1_db"
    with pytest.raises(checkm1.CheckM1Error, match="overlap"):
        checkm1.validate_destination(
            destination, [tmp_path / "results"], tmp_path / "environment"
        )
    actual_environment = tmp_path / "environment"
    actual_environment.mkdir()
    environment_alias = tmp_path / "environment-alias"
    environment_alias.symlink_to(actual_environment, target_is_directory=True)
    first = checkm1.environment_lock_path(actual_environment)
    second = checkm1.environment_lock_path(environment_alias)
    assert first == second == actual_environment / ".meew-checkm1-environment.lock"
    assert "/tmp/meew-checkm1-environment" not in str(first)
    assert str(os.getuid()) not in first.name


@pytest.mark.parametrize("lock_kind", ["cache", "environment"])
def test_cache_and_environment_locks_serialize_local_processes(tmp_path, lock_kind):
    if lock_kind == "cache":
        lock = checkm1.cache_lock_path(tmp_path / "database")
    else:
        environment = tmp_path / "environment"
        environment.mkdir()
        lock = checkm1.environment_lock_path(environment)
    events = tmp_path / (lock_kind + ".events")
    script = """
import importlib.util
import pathlib
import sys
import time
spec = importlib.util.spec_from_file_location("checkm1_lock_worker", sys.argv[1])
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
with module._exclusive_lock(pathlib.Path(sys.argv[2]), sys.argv[3]):
    with open(sys.argv[4], "a", encoding="utf-8") as handle:
        handle.write(sys.argv[5] + "-start\\n")
    time.sleep(0.2)
    with open(sys.argv[4], "a", encoding="utf-8") as handle:
        handle.write(sys.argv[5] + "-end\\n")
"""
    environment = dict(os.environ)
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    workers = []
    for label in ("one", "two"):
        workers.append(
            subprocess.Popen(
                [
                    sys.executable,
                    "-c",
                    script,
                    str(HELPER),
                    str(lock),
                    lock_kind,
                    str(events),
                    label,
                ],
                env=environment,
            )
        )
        time.sleep(0.05)
    assert [worker.wait() for worker in workers] == [0, 0]
    lines = events.read_text(encoding="utf-8").splitlines()
    assert lines in (
        ["one-start", "one-end", "two-start", "two-end"],
        ["two-start", "two-end", "one-start", "one-end"],
    )


def test_cancellation_before_launch_never_executes_command(tmp_path):
    marker = tmp_path / "executed"
    with checkm1._SignalForwarder() as forwarder:
        forwarder._handle(signal.SIGTERM, None)
        returncode = forwarder.run(
            [sys.executable, "-c", "open(%r, 'w').write('bad')" % str(marker)]
        )
    assert returncode == 143
    assert not marker.exists()


def test_lock_blocked_cancellation_aborts_promptly_and_releases_waiter(tmp_path):
    lock = tmp_path / "blocked.lock"
    ready = tmp_path / "holder.ready"
    holder_code = (
        "import fcntl,sys,time; "
        "h=open(sys.argv[1],'a+'); fcntl.flock(h.fileno(),fcntl.LOCK_EX); "
        "open(sys.argv[2],'w').write('ready'); time.sleep(0.8)"
    )
    holder = subprocess.Popen(
        [sys.executable, "-c", holder_code, str(lock), str(ready)]
    )
    deadline = time.monotonic() + 2
    while not ready.exists() and time.monotonic() < deadline:
        time.sleep(0.01)
    assert ready.exists()

    with checkm1._SignalForwarder() as forwarder:
        timer = threading.Timer(0.15, lambda: os.kill(os.getpid(), signal.SIGTERM))
        timer.start()
        started = time.monotonic()
        with pytest.raises(checkm1.CheckM1Cancelled):
            with checkm1._exclusive_lock(lock, "cache", cancellation=forwarder):
                raise AssertionError("blocked waiter acquired the lock")
        timer.join()
        assert time.monotonic() - started < 0.7
    assert holder.wait() == 0
    with checkm1._exclusive_lock(lock, "cache"):
        pass


def test_active_child_cancellation_stops_group_no_continuation_and_releases_lock(
    tmp_path,
):
    lock = tmp_path / "active.lock"
    started = tmp_path / "started"
    continued = tmp_path / "continued"
    child_code = (
        "import pathlib,sys,time; "
        "pathlib.Path(sys.argv[1]).write_text('started'); "
        "time.sleep(5); pathlib.Path(sys.argv[2]).write_text('continued')"
    )
    with checkm1._SignalForwarder() as forwarder:
        with checkm1._exclusive_lock(lock, "environment", cancellation=forwarder):
            timer = threading.Timer(0.2, lambda: os.kill(os.getpid(), signal.SIGTERM))
            timer.start()
            returncode = forwarder.run(
                [sys.executable, "-c", child_code, str(started), str(continued)]
            )
            timer.join()
    assert returncode == 143
    assert started.exists()
    assert not continued.exists()
    with checkm1._exclusive_lock(lock, "environment"):
        pass


@pytest.mark.skipif(
    not hasattr(signal, "pthread_sigmask"), reason="POSIX signal masks unavailable"
)
def test_prompt_group_leader_keeps_descendant_protected_with_usable_signal_mask(
    tmp_path,
):
    lock = tmp_path / "descendant.lock"
    ready = tmp_path / "descendant.ready"
    terminated = tmp_path / "descendant.terminated"
    continued = tmp_path / "descendant.continued"
    mask = tmp_path / "descendant.mask"
    descendant_code = (
        "import pathlib,signal,sys,time; "
        "ready=pathlib.Path(sys.argv[1]); terminated=pathlib.Path(sys.argv[2]); "
        "continued=pathlib.Path(sys.argv[3]); mask=pathlib.Path(sys.argv[4]); "
        "blocked=signal.pthread_sigmask(signal.SIG_BLOCK, []); "
        "mask.write_text(str(signal.SIGTERM in blocked)); "
        "signal.signal(signal.SIGTERM, lambda signum, frame: "
        "(terminated.write_text('terminated'), sys.exit(0))); "
        "ready.write_text('ready'); time.sleep(5); continued.write_text('continued')"
    )
    leader_code = (
        "import subprocess,sys; "
        "subprocess.Popen([sys.executable,'-c',sys.argv[1],*sys.argv[2:]], "
        "close_fds=True)"
    )

    def cancel_when_ready():
        deadline = time.monotonic() + 2
        while not ready.exists() and time.monotonic() < deadline:
            time.sleep(0.01)
        os.kill(os.getpid(), signal.SIGTERM)

    with checkm1._SignalForwarder() as forwarder:
        with checkm1._exclusive_lock(lock, "environment", cancellation=forwarder):
            canceller = threading.Thread(target=cancel_when_ready)
            canceller.start()
            returncode = forwarder.run(
                [
                    sys.executable,
                    "-c",
                    leader_code,
                    descendant_code,
                    str(ready),
                    str(terminated),
                    str(continued),
                    str(mask),
                ]
            )
            canceller.join()
            assert terminated.exists()
            assert not continued.exists()
        with checkm1._exclusive_lock(lock, "environment"):
            pass

    assert returncode == 143
    assert mask.read_text(encoding="utf-8") == "False"


def test_registration_hardlink_isolation_preserves_other_environment(tmp_path):
    first = tmp_path / "env-one" / "checkm" / "DATA_CONFIG"
    second = tmp_path / "env-two" / "checkm" / "DATA_CONFIG"
    first.parent.mkdir(parents=True)
    second.parent.mkdir(parents=True)
    snapshot = b'{"dataRoot": "/prior"}'
    first.write_bytes(snapshot)
    os.link(str(first), str(second))
    assert first.stat().st_ino == second.stat().st_ino

    assert checkm1._isolate_registration_hardlinks(first, snapshot) is True
    first.write_text('{"dataRoot": "/selected"}', encoding="utf-8")
    assert second.read_bytes() == snapshot
    assert first.stat().st_ino != second.stat().st_ino
    assert first.stat().st_nlink == 1


def test_registration_updates_only_data_config_not_validated_database(tmp_path):
    database = _ready_database(tmp_path / "database")
    manifest = database / ".dmanifest"
    manifest_before = (manifest.read_bytes(), manifest.stat().st_mtime_ns)
    config_path = tmp_path / "environment" / "checkm" / "DATA_CONFIG"
    config_path.parent.mkdir(parents=True)
    config_path.write_text(
        json.dumps(
            {
                "dataRoot": "/prior",
                "remoteManifestURL": "https://example.invalid/",
                "manifestType": "CheckM",
                "localManifestName": ".dmanifest",
                "remoteManifestName": ".dmanifest",
            }
        ),
        encoding="utf-8",
    )
    marker_root = _marker_root(tmp_path / "markers")

    class RegistrationForwarder(_Forwarder):
        def run_capture(self, command, environ=None):
            del command, environ
            return 0, json.dumps({"root": str(database), "counts": [104, 149]}), ""

    checkm1._register_database(
        database,
        config_path,
        marker_root,
        RegistrationForwarder(),
        {},
    )
    assert json.loads(config_path.read_text(encoding="utf-8"))["dataRoot"] == str(
        database
    )
    assert (manifest.read_bytes(), manifest.stat().st_mtime_ns) == manifest_before


def test_registration_failure_restores_config_and_preserves_valid_cache(
    tmp_path, monkeypatch
):
    destination = _ready_database(tmp_path / "database")
    prefix = tmp_path / "environment"
    config_path = (
        prefix / "lib" / "python3.7" / "site-packages" / "checkm" / "DATA_CONFIG"
    )
    config_path.parent.mkdir(parents=True)
    snapshot = b'{"dataRoot": "/prior"}'
    config_path.write_bytes(snapshot)
    marker_root = _marker_root(tmp_path / "markers")

    monkeypatch.setattr(
        checkm1,
        "validate_native_installation",
        lambda *args, **kwargs: {
            "hook": tmp_path / "hook",
            "marker_root": marker_root,
        },
    )
    monkeypatch.setattr(checkm1, "_checkm_data_config_path", lambda: config_path)
    monkeypatch.setattr(checkm1, "_prepare_database", lambda *args: None)

    def corrupt_then_fail(database, path, *args):
        del database, args
        path.write_text("corrupt", encoding="utf-8")
        raise checkm1.CheckM1Error("registration denied")

    monkeypatch.setattr(checkm1, "_register_database", corrupt_then_fail)
    with pytest.raises(checkm1.CheckM1Error, match="registration denied"):
        checkm1.prepare_and_run(
            destination,
            True,
            [],
            ["never-run"],
            environment_prefix=prefix,
        )
    assert config_path.read_bytes() == snapshot
    checkm1.validate_database(destination)


def test_corrupt_registration_and_restore_permissions_fail_before_command(
    tmp_path, monkeypatch
):
    config_path = tmp_path / "environment" / "checkm" / "DATA_CONFIG"
    config_path.parent.mkdir(parents=True)
    config_path.write_text("not json", encoding="utf-8")
    with pytest.raises(checkm1.CheckM1Error, match="corrupt"):
        checkm1._registration_snapshot(config_path)

    config_path.write_text('{"dataRoot": "/prior"}', encoding="utf-8")

    def denied(*args, **kwargs):
        del args, kwargs
        raise PermissionError("read-only environment")

    monkeypatch.setattr(checkm1.tempfile, "mkstemp", denied)
    with pytest.raises(checkm1.CheckM1Error, match="restore prior.*read-only"):
        checkm1._restore_registration(config_path, config_path.read_bytes())


def test_successful_preflight_smoke_launches_only_configured_cheap_command(
    tmp_path, monkeypatch
):
    destination = _ready_database(tmp_path / "database")
    prefix = tmp_path / "environment"
    config_path = prefix / "checkm" / "DATA_CONFIG"
    config_path.parent.mkdir(parents=True)
    config_path.write_text('{"dataRoot": "/prior"}', encoding="utf-8")
    marker_root = _marker_root(tmp_path / "markers")
    launched = tmp_path / "cheap-command-launched"

    monkeypatch.setattr(
        checkm1,
        "validate_native_installation",
        lambda *args, **kwargs: {
            "hook": tmp_path / "hook",
            "marker_root": marker_root,
        },
    )
    monkeypatch.setattr(checkm1, "_checkm_data_config_path", lambda: config_path)
    monkeypatch.setattr(checkm1, "_prepare_database", lambda *args: None)

    def register(database, path, *args):
        del args
        path.write_text(json.dumps({"dataRoot": str(database)}), encoding="utf-8")

    monkeypatch.setattr(checkm1, "_register_database", register)
    result = checkm1.prepare_and_run(
        destination,
        True,
        [],
        [
            sys.executable,
            "-c",
            "from pathlib import Path; Path(%r).write_text('cheap')" % str(launched),
        ],
        environment_prefix=prefix,
    )
    assert result == 0
    assert launched.read_text(encoding="utf-8") == "cheap"


def test_rule_declares_shared_sources_for_both_binners_and_wraps_cleanup():
    source = RULES.read_text(encoding="utf-8")
    comebin = source.split("\n\nrule metabinner_depth:", 1)[0]
    metabinner = source.split("\n\nrule metabinner:\n", 1)[1].split(
        "\n\nrule metabinner_bins:", 1
    )[0]
    for rule in (comebin, metabinner):
        assert "checkm1_runner=CHECKM1_RUNNER" in rule
        assert "checkm1_helper=CHECKM1_HELPER" in rule
    assert 'atomic_script=workflow.source_path("../scripts/run_comebin_atomic.sh")' in (
        comebin
    )
    assert "bash {input.atomic_script:q}" in comebin
    assert "--checkm-runner {input.checkm1_runner:q}" in comebin
    assert "--checkm-helper {input.checkm1_helper:q}" in comebin
    assert comebin.index("--database {params.checkm1_database:q}") < comebin.index(
        "-- {params.extra}"
    )
    assert "python {input.checkm1_runner:q}" in metabinner
    assert "--helper-source {input.checkm1_helper:q}" in metabinner
    assert "--tool metabinner" in metabinner
    assert '--tool-install-root "$install_root"' in metabinner
    assert '--cleanup-path "$output_root"' in metabinner
    assert 'rm -rf "$output_root"' not in metabinner
    assert 'ambient_install="$(printenv METABINNER_HOME || true)"' in metabinner
    native_arguments = [
        '-- bash "$runner"',
        '-a "$contigs"',
        '-o "$output_root"',
        '-d "$coverage"',
        '-k "$composition"',
        '-p "$install_root"',
        "-t {threads}",
        "-s {params.scale:q}",
    ]
    assert [metabinner.index(argument) for argument in native_arguments] == sorted(
        metabinner.index(argument) for argument in native_arguments
    )
    assert RUNNER.is_file()
    assert COMEBIN_ATOMIC.is_file()
    assert METABINNER_PROFILES.is_file()


def test_standalone_dry_run_quotes_sources_without_setup_side_effects(tmp_path):
    workdir = tmp_path / "workflow with spaces"
    workdir.mkdir()
    inputs = workdir / "input data"
    inputs.mkdir()
    contigs = inputs / "assembly file.fa"
    bam = inputs / "alignment file.bam"
    contigs.write_text(">contig\nACGT\n", encoding="utf-8")
    bam.write_text("fixture\n", encoding="utf-8")
    samples = workdir / "samples.tsv"
    samples.write_text(
        "sample\tcontigs\tbam\n" "space_sample\t%s\t%s\n" % (contigs, bam),
        encoding="utf-8",
    )
    resource_root = workdir / "resources not created"
    config = yaml.safe_load(
        (ROOT / "config" / "config.yaml").read_text(encoding="utf-8")
    )
    config["sample_sheet"] = str(samples)
    config["resources_root"] = str(resource_root)
    config_path = workdir / "config.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")

    result = subprocess.run(
        [
            "snakemake",
            "--snakefile",
            str(ROOT / "workflow" / "Snakefile"),
            "--directory",
            str(workdir),
            "--configfile",
            str(config_path),
            "--replace-workflow-config",
            "--cores",
            "1",
            "--dry-run",
            "--printshellcmds",
        ],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    rendered = result.stdout + result.stderr
    assert str(resource_root / "checkm1_db") in rendered
    assert "run_comebin_with_checkm.py" in rendered
    assert "run_comebin_atomic.sh" in rendered
    assert "source-cache" in rendered
    assert not resource_root.exists()


def test_standalone_four_binner_refinement_dag_keeps_metabinner_sources():
    result = subprocess.run(
        [
            "snakemake",
            "--snakefile",
            str(ROOT / "workflow" / "Snakefile"),
            "--directory",
            str(ROOT),
            "--configfile",
            str(ROOT / ".test" / "config" / "config.refinement.yaml"),
            "--replace-workflow-config",
            "--cores",
            "1",
            "--dry-run",
            "--printshellcmds",
        ],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    rendered = result.stdout + result.stderr
    for rule in (
        "comebin",
        "semibin2",
        "metabat2",
        "metabinner_depth",
        "metabinner_profiles",
        "metabinner",
        "metabinner_bins",
    ):
        assert "rule %s:" % rule in rendered
    assert "refine_das_tool" in rendered
    metabinner_job = rendered.split("rule metabinner:", 1)[1].split("\n\nrule ", 1)[0]
    assert "run_comebin_with_checkm.py" in metabinner_job
    assert "common/checkm1.py" in metabinner_job
    assert "--tool metabinner" in metabinner_job
    profiles_job = rendered.split("rule metabinner_profiles:", 1)[1].split(
        "\n\nrule ", 1
    )[0]
    assert "prepare_metabinner_profiles.sh" in profiles_job
    assert "source-cache" in profiles_job


def test_thin_runner_preserves_arguments_with_spaces(tmp_path):
    helper_dir = tmp_path / "source cache with spaces"
    helper_dir.mkdir()
    fake_helper = helper_dir / "fake helper.py"
    fake_helper.write_text(
        "import json\n"
        "def main(arguments=None):\n"
        "    print(json.dumps(arguments))\n"
        "    return 0\n",
        encoding="utf-8",
    )
    result = subprocess.run(
        [
            sys.executable,
            str(RUNNER),
            "--helper-source",
            str(fake_helper),
            "--database",
            str(tmp_path / "database with spaces"),
            "--",
            "command with spaces",
            "argument with spaces",
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert json.loads(result.stdout) == [
        "--database",
        str(tmp_path / "database with spaces"),
        "--",
        "command with spaces",
        "argument with spaces",
    ]
