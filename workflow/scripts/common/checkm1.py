"""Prepare CheckM 1 data safely and run one CheckM-consuming binner.

This module is intentionally compatible with the Python 3.7 runtimes shipped
by the pinned COMEBin and MetaBinner environments.  Importing it has no setup,
registration, lock, or download side effects.
"""

from __future__ import print_function

import argparse
import contextlib
import errno
import fcntl
import hashlib
import importlib.util
import json
import os
import re
import shutil
import signal
import stat
import subprocess
import sys
import tempfile
import time
from pathlib import Path

OFFICIAL_ARCHIVE_URL = (
    "https://data.ace.uq.edu.au/public/CheckM_databases/"
    "checkm_data_2015_01_16.tar.gz"
)
OFFICIAL_ARCHIVE_SHA256 = (
    "971ec469348bd6c3d9eb96142f567f12443310fa06c1892643940f35f86ac92c"
)
REQUIRED_DATABASE_FILES = (
    ".dmanifest",
    "hmms/checkm.hmm",
    "hmms/checkm.hmm.ssi",
    "pfam/Pfam-A.hmm.dat",
)
REQUIRED_SCORING_EXECUTABLES = (
    "checkm",
    "hmmsearch",
    "hmmfetch",
    "hmmalign",
    "prodigal",
)
REQUIRED_SETUP_EXECUTABLES = ("bash", "sha256sum", "tar", "wget")
MARKER_FILES = ("bacteria.ms", "archaea.ms")
SUPPORTED_TOOLS = ("comebin", "metabinner")
MAX_MANIFEST_BYTES = 1024 * 1024
MAX_FULL_HASH_BYTES = 8 * 1024 * 1024
PROBE_BYTES = 64 * 1024
LOCK_POLL_SECONDS = 0.1
CHILD_TERMINATION_GRACE_SECONDS = 3.0
_SIGNALS = (signal.SIGHUP, signal.SIGINT, signal.SIGTERM)


class CheckM1Error(RuntimeError):
    """A safe CheckM 1 setup or registration operation failed."""


class CheckM1Cancelled(CheckM1Error):
    """The wrapper received a termination signal and aborted safely."""

    def __init__(self, signum):
        self.signum = signum
        CheckM1Error.__init__(self, "CheckM 1 runtime cancelled by signal %s." % signum)

    @property
    def returncode(self):
        return 128 + self.signum


class CheckM1RegistrationProtectionError(CheckM1Error):
    """Registration rollback failed, so staging must remain available."""

    preserve_staging = True


def expand_path(value, environ=None):
    """Expand variables and ``~`` while preserving relative path semantics."""
    if value is None or not str(value).strip():
        return ""
    environment = os.environ if environ is None else environ
    expanded = (
        os.path.expandvars(str(value))
        if environment is os.environ
        else _expandvars(str(value), environment)
    )
    home = environment.get("HOME", "")
    if expanded == "~" and home:
        return home
    if expanded.startswith("~/") and home:
        return os.path.join(home, expanded[2:])
    return os.path.expanduser(expanded)


def _expandvars(value, environ):
    """Expand common POSIX variable forms against an explicit environment."""
    result = value
    for name, configured in environ.items():
        result = result.replace("${%s}" % name, configured)
        result = result.replace("$%s" % name, configured)
    return result


def resolve_checkm1_database(explicit="", module_root="", environ=None):
    """Resolve the exact CheckM 1 database directory without creating it."""
    environment = os.environ if environ is None else environ
    if explicit is not None and str(explicit).strip():
        return expand_path(explicit, environment)
    if module_root is not None and str(module_root).strip():
        root = Path(expand_path(module_root, environment))
    elif environment.get("MEEW_RESOURCES", "").strip():
        root = Path(expand_path(environment["MEEW_RESOURCES"], environment))
    elif environment.get("XDG_DATA_HOME", "").strip():
        root = Path(expand_path(environment["XDG_DATA_HOME"], environment)) / "meew"
    else:
        root = Path(expand_path("~/resources", environment))
    return str(root / "checkm1_db")


def _nonempty_readable_file(path):
    try:
        if not path.is_file() or path.stat().st_size <= 0:
            return False
        with path.open("rb") as handle:
            return bool(handle.read(1))
    except OSError:
        return False


def _sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _bounded_edges(path):
    size = path.stat().st_size
    with path.open("rb") as handle:
        first = handle.read(min(PROBE_BYTES, size))
        handle.seek(max(0, size - PROBE_BYTES))
        last = handle.read(PROBE_BYTES)
    return first, last


def _manifest_entries(manifest):
    try:
        if manifest.stat().st_size > MAX_MANIFEST_BYTES:
            raise ValueError("manifest exceeds the bounded parser limit")
        lines = manifest.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeError) as error:
        raise ValueError("manifest is unreadable: %s" % error)
    if not lines or lines[0] != (
        "##CheckM##\tData manifest created by ScreamingBackpack version 0.2.3"
    ):
        raise ValueError("manifest header is not the pinned CheckM 1 format")

    entries = {}
    for line_number, line in enumerate(lines[1:], 2):
        fields = line.split("\t")
        if len(fields) != 3:
            raise ValueError("malformed manifest row %d" % line_number)
        relative, digest, size_text = fields
        try:
            size = int(size_text)
        except ValueError:
            raise ValueError("invalid manifest size on row %d" % line_number)
        if size < 0 or (digest != "-" and not re.match(r"^[0-9a-f]{64}$", digest)):
            raise ValueError("invalid manifest metadata on row %d" % line_number)
        value = (digest, size)
        if relative in entries and entries[relative] != value:
            raise ValueError("conflicting duplicate manifest row for %s" % relative)
        entries[relative] = value
    return entries


def _format_problems(root, entries):
    problems = []
    model = root / "hmms" / "checkm.hmm"
    index = root / "hmms" / "checkm.hmm.ssi"
    pfam = root / "pfam" / "Pfam-A.hmm.dat"

    for relative in REQUIRED_DATABASE_FILES[1:]:
        path = root / relative
        if relative not in entries:
            problems.append("%s is absent from .dmanifest" % relative)
            continue
        if not _nonempty_readable_file(path):
            problems.append("%s is missing, unreadable, or empty" % relative)
            continue
        digest, expected_size = entries[relative]
        try:
            actual_size = path.stat().st_size
        except OSError as error:
            problems.append("%s cannot be inspected: %s" % (relative, error))
            continue
        if actual_size != expected_size:
            problems.append(
                "%s size is %d but .dmanifest requires %d"
                % (relative, actual_size, expected_size)
            )
            continue
        if actual_size <= MAX_FULL_HASH_BYTES:
            try:
                actual_digest = _sha256(path)
            except OSError as error:
                problems.append("%s cannot be hashed: %s" % (relative, error))
                continue
            if actual_digest != digest:
                problems.append("%s SHA-256 does not match .dmanifest" % relative)

    if _nonempty_readable_file(model):
        try:
            first, last = _bounded_edges(model)
            required = (b"HMMER3/", b"\nNAME ", b"\nACC ", b"\nLENG ", b"\nHMM ")
            if not first.startswith(required[0]) or any(
                token not in first for token in required[1:]
            ):
                problems.append("hmms/checkm.hmm has an invalid HMMER model header")
            if not last.rstrip().endswith(b"//"):
                problems.append("hmms/checkm.hmm has an invalid/truncated HMMER tail")
        except OSError as error:
            problems.append("hmms/checkm.hmm format cannot be inspected: %s" % error)

    if _nonempty_readable_file(index):
        try:
            first, _ = _bounded_edges(index)
            if (
                not first.startswith(b"\xd3\xd3\xc9\xb3")
                or b"checkm.hmm\x00" not in first
            ):
                problems.append("hmms/checkm.hmm.ssi is not a CheckM HMM SSI index")
        except OSError as error:
            problems.append(
                "hmms/checkm.hmm.ssi format cannot be inspected: %s" % error
            )

    if _nonempty_readable_file(pfam):
        try:
            first, last = _bounded_edges(pfam)
            if not first.startswith(b"# STOCKHOLM 1.0\n"):
                problems.append("pfam/Pfam-A.hmm.dat lacks a Stockholm header")
            if b"#=GF AC   PF" not in first or b"#=GF GA   " not in first:
                problems.append(
                    "pfam/Pfam-A.hmm.dat lacks PFAM accession/threshold data"
                )
            if not last.rstrip().endswith(b"//"):
                problems.append("pfam/Pfam-A.hmm.dat has an invalid/truncated tail")
        except OSError as error:
            problems.append(
                "pfam/Pfam-A.hmm.dat format cannot be inspected: %s" % error
            )
    return problems


def database_problems(database):
    """Return bounded, meaningful readiness defects for a CheckM 1 cache."""
    root = Path(database)
    manifest = root / ".dmanifest"
    if not _nonempty_readable_file(manifest):
        return [".dmanifest is missing, unreadable, or empty"]
    try:
        entries = _manifest_entries(manifest)
    except ValueError as error:
        return [".dmanifest is invalid: %s" % error]
    return _format_problems(root, entries)


def validate_database(database):
    """Require format-consistent native CheckM analyze/qa payloads."""
    problems = database_problems(database)
    if problems:
        raise CheckM1Error(
            "CheckM 1 database is incomplete or corrupt at %s: %s"
            % (database, "; ".join(problems))
        )


def _which(executable, environ):
    return shutil.which(executable, path=environ.get("PATH"))


def _command_path(value, environ):
    resolved = _which(value, environ)
    if resolved is not None:
        return Path(resolved).resolve()
    candidate = Path(value)
    if candidate.is_file():
        return candidate.resolve()
    return None


def _marker_file_problems(path, expected_domain):
    if not _nonempty_readable_file(path):
        return ["%s is missing, unreadable, or empty" % path]
    try:
        if path.stat().st_size > MAX_MANIFEST_BYTES:
            return ["%s exceeds the bounded marker parser limit" % path]
        lines = [
            line
            for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
    except (OSError, UnicodeError) as error:
        return ["%s is unreadable: %s" % (path, error)]
    if len(lines) != 2 or lines[0] != "# [Taxon Marker File]":
        return ["%s is not a two-row CheckM taxon marker file" % path]
    fields = lines[1].split("\t")
    if len(fields) != 6 or fields[0] != expected_domain or fields[3] != expected_domain:
        return ["%s has malformed %s marker metadata" % (path, expected_domain)]
    try:
        if int(fields[1]) < 1 or int(fields[4]) < 1:
            raise ValueError("nonpositive marker metadata")
        int(fields[2])
    except ValueError:
        return ["%s has nonnumeric marker metadata" % path]
    expression = fields[5]
    marker_ids = set(re.findall(r"(?:PF\d{5}(?:\.\d+)?|TIGR\d{5})", expression))
    if not expression.startswith("[set([") or not expression.endswith(")]"):
        return ["%s has malformed marker-set structure" % path]
    if len(marker_ids) < 40:
        return ["%s has too few recognizable PFAM/TIGRFAM markers" % path]
    return []


def _tool_installation(tool, command, environ, tool_install_root=""):
    if tool not in SUPPORTED_TOOLS:
        raise CheckM1Error("Unsupported CheckM-consuming tool: %s" % tool)
    if not command:
        raise CheckM1Error("No %s command was supplied." % tool)

    if tool == "comebin":
        runner = _command_path(command[0], environ)
        if runner is None:
            raise CheckM1Error("COMEBin runner is unavailable: %s" % command[0])
        marker_root = runner.parent / "COMEBin" / "scripts" / "checkm_ms"
        install_root = runner.parent / "COMEBin"
    else:
        if not tool_install_root or not str(tool_install_root).strip():
            raise CheckM1Error(
                "MetaBinner requires its resolved tool installation root."
            )
        install_root = _canonical(tool_install_root)
        runner = install_root / "run_metabinner.sh"
        command_runner = None
        first = _command_path(command[0], environ)
        if first is not None and first.name == "bash" and len(command) > 1:
            command_runner = _command_path(command[1], environ)
        elif first is not None:
            command_runner = first
        if not runner.is_file() or command_runner is None:
            raise CheckM1Error(
                "MetaBinner runner is unavailable below the resolved tool root: %s"
                % runner
            )
        try:
            same_runner = os.path.samefile(str(runner), str(command_runner))
        except OSError:
            same_runner = False
        if not same_runner:
            raise CheckM1Error(
                "MetaBinner command runner %s is not %s; the shell launcher must "
                "not be treated as the tool installation root."
                % (command_runner, runner)
            )
        if (
            not (install_root / "scripts").is_dir()
            or not (install_root / "auxiliary").is_dir()
        ):
            raise CheckM1Error(
                "Incomplete MetaBinner installation at %s" % install_root
            )
        marker_root = install_root / "scripts" / "checkm_ms"

    marker_problems = []
    for marker, domain in zip(MARKER_FILES, ("Bacteria", "Archaea")):
        marker_problems.extend(_marker_file_problems(marker_root / marker, domain))
    if marker_problems:
        raise CheckM1Error(
            "%s bundled marker readiness failed: %s"
            % (tool, "; ".join(marker_problems))
        )
    return {
        "tool": tool,
        "runner": runner,
        "install_root": install_root,
        "marker_root": marker_root,
    }


def validate_native_installation(
    tool, command, environ=None, environment_prefix=None, tool_install_root=""
):
    """Validate the installed hook, scoring tools, runner, and marker sets."""
    environment = dict(os.environ if environ is None else environ)
    prefix = _canonical(environment_prefix or sys.prefix)
    hook = prefix / "bin" / ".checkm-genome-post-link.sh"
    if not _nonempty_readable_file(hook):
        raise CheckM1Error(
            "Installed CheckM Bioconda post-link hook is missing or unreadable: %s"
            % hook
        )

    hook_text = hook.read_text(encoding="utf-8")
    required_hook_text = (
        "CHECKM_DATA_DIR",
        "CHECKM_DATA_PATH",
        OFFICIAL_ARCHIVE_URL,
        OFFICIAL_ARCHIVE_SHA256,
        "sha256sum -c",
        "checkm data setRoot",
    )
    missing_hook_contract = [
        token for token in required_hook_text if token not in hook_text
    ]
    if missing_hook_contract:
        raise CheckM1Error(
            "Installed CheckM post-link hook does not match the pinned native setup "
            "contract; missing token(s): %s" % ", ".join(missing_hook_contract)
        )

    missing_executables = [
        executable
        for executable in REQUIRED_SETUP_EXECUTABLES + REQUIRED_SCORING_EXECUTABLES
        if not _which(executable, environment)
    ]
    if missing_executables:
        raise CheckM1Error(
            "Required CheckM setup/scoring executable(s) are unavailable on PATH: %s"
            % ", ".join(missing_executables)
        )

    installation = _tool_installation(
        tool, command, environment, tool_install_root=tool_install_root
    )
    installation["hook"] = hook
    return installation


def _canonical(path):
    return Path(path).expanduser().absolute().resolve(strict=False)


def paths_overlap(first, second):
    """Return whether either canonical path contains the other."""
    first_path = _canonical(first)
    second_path = _canonical(second)
    return (
        first_path == second_path
        or first_path in second_path.parents
        or second_path in first_path.parents
    )


def _contains_whitespace(path):
    return bool(re.search(r"\s", str(path)))


def validate_native_paths(database, tool, tool_install_root="", native_paths=None):
    """Reject paths unsupported by native CheckM/MetaBinner shell commands."""
    candidates = [("CheckM 1 database", database)]
    if tool == "metabinner":
        candidates.append(("MetaBinner installation", tool_install_root))
        candidates.extend(
            ("MetaBinner native input/output", value)
            for value in (native_paths or [])
            if value
        )
    unsupported = [
        "%s %s" % (label, value)
        for label, value in candidates
        if value and _contains_whitespace(value)
    ]
    if unsupported:
        raise CheckM1Error(
            "Whitespace-containing native paths are unsupported because CheckM 1 "
            "HMM commands and MetaBinner 1.4.4 shell code do not quote every path: %s"
            % "; ".join(unsupported)
        )


def validate_destination(database, protected_paths, environment_prefix=None):
    """Reject destinations that could consume or be consumed by protected data."""
    requested_destination = Path(database).expanduser().absolute()
    if requested_destination.is_symlink():
        raise CheckM1Error(
            "CheckM 1 database destination must not be a symlink: %s"
            % requested_destination
        )
    destination = requested_destination.resolve(strict=False)
    if destination == destination.parent:
        raise CheckM1Error(
            "CheckM 1 database destination cannot be a filesystem root: %s"
            % destination
        )
    protected = [("active environment", environment_prefix or sys.prefix)]
    protected.extend(
        ("workflow input/result/workspace", value) for value in protected_paths if value
    )
    overlaps = [
        "%s %s" % (label, value)
        for label, value in protected
        if paths_overlap(destination, value)
    ]
    if overlaps:
        raise CheckM1Error(
            "Refusing destructive CheckM 1 destination overlap at %s with: %s"
            % (destination, "; ".join(overlaps))
        )
    return destination


@contextlib.contextmanager
def _exclusive_lock(path, description, cancellation=None):
    path = Path(path)
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        handle = path.open("a+")
        if not stat.S_ISREG(os.fstat(handle.fileno()).st_mode):
            raise OSError("lock is not a regular file")
    except OSError as error:
        raise CheckM1Error(
            "Cannot establish shared %s lock %s: %s. Ensure the shared "
            "environment/cache location is writable and visible to every worker."
            % (description, path, error)
        )
    acquired = False
    try:
        while not acquired:
            if cancellation is not None:
                cancellation.raise_if_cancelled()
            try:
                fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
                acquired = True
            except OSError as error:
                if error.errno not in (errno.EACCES, errno.EAGAIN):
                    raise CheckM1Error(
                        "Cannot acquire shared %s lock %s: %s"
                        % (description, path, error)
                    )
                time.sleep(LOCK_POLL_SECONDS)
        if cancellation is not None:
            cancellation.raise_if_cancelled()
        yield
    finally:
        try:
            if acquired:
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
        finally:
            handle.close()


def environment_lock_path(environment_prefix=None):
    """Return a shared-filesystem lock path tied to one canonical environment."""
    prefix = _canonical(environment_prefix or sys.prefix)
    return prefix / ".meew-checkm1-environment.lock"


def cache_lock_path(database):
    destination = _canonical(database)
    return destination.parent / (".%s.meew-checkm1.lock" % destination.name)


def _checkm_data_config_path():
    spec = importlib.util.find_spec("checkm")
    if spec is None or not spec.submodule_search_locations:
        raise CheckM1Error("The installed CheckM Python package cannot be located.")
    return Path(next(iter(spec.submodule_search_locations))) / "DATA_CONFIG"


def _registration_snapshot(config_path):
    try:
        if config_path.is_symlink() or not config_path.is_file():
            raise ValueError("DATA_CONFIG must be a regular, non-symlink file")
        payload = config_path.read_bytes()
        values = json.loads(payload.decode("utf-8"))
        if not isinstance(values.get("dataRoot"), str):
            raise ValueError("dataRoot is not a string")
        return payload
    except (OSError, UnicodeError, ValueError, TypeError) as error:
        raise CheckM1Error(
            "CheckM DATA_CONFIG is unreadable or corrupt at %s: %s"
            % (config_path, error)
        )


def _validate_registration_location(config_path, environment_prefix):
    prefix = _canonical(environment_prefix)
    requested = Path(config_path).expanduser().absolute()
    if requested.is_symlink():
        raise CheckM1Error("CheckM DATA_CONFIG must not be a symlink: %s" % requested)
    resolved = requested.resolve(strict=True)
    if prefix != resolved and prefix not in resolved.parents:
        raise CheckM1Error(
            "CheckM DATA_CONFIG %s is outside the protected environment %s. "
            "Refusing cross-environment registration." % (resolved, prefix)
        )


def _atomic_registration_write(config_path, payload, operation):
    temporary = None
    try:
        stat_result = config_path.stat()
        descriptor, temporary_name = tempfile.mkstemp(
            prefix=".%s.%s-" % (config_path.name, operation),
            dir=str(config_path.parent),
        )
        temporary = Path(temporary_name)
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.chmod(str(temporary), stat_result.st_mode & 0o777)
        os.replace(str(temporary), str(config_path))
    except OSError as error:
        raise CheckM1Error(
            "Failed to %s CheckM DATA_CONFIG at %s: %s"
            % (operation, config_path, error)
        )
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def _restore_registration(config_path, payload):
    _atomic_registration_write(config_path, payload, "restore prior")


def _isolate_registration_hardlinks(config_path, snapshot):
    """Detach DATA_CONFIG before native code can mutate shared hard links."""
    try:
        links = config_path.stat().st_nlink
    except OSError as error:
        raise CheckM1Error("Cannot inspect CheckM DATA_CONFIG links: %s" % error)
    if links <= 1:
        return False
    _atomic_registration_write(config_path, snapshot, "isolate shared")
    try:
        if config_path.stat().st_nlink != 1 or config_path.read_bytes() != snapshot:
            raise CheckM1Error(
                "CheckM DATA_CONFIG hard-link isolation could not be verified at %s"
                % config_path
            )
    except OSError as error:
        raise CheckM1Error(
            "CheckM DATA_CONFIG hard-link isolation verification failed: %s" % error
        )
    print(
        "INFO: detached CheckM DATA_CONFIG from %d shared hard links before registration"
        % links
    )
    return True


class _SignalForwarder(object):
    def __init__(self):
        self.process = None
        self.process_group = None
        self.received = None
        self.previous = {}

    def _process_group_exists(self):
        if self.process_group is None:
            return False
        try:
            os.killpg(self.process_group, 0)
        except ProcessLookupError:
            return False
        except PermissionError:
            return True
        return True

    def _forward(self, signum):
        if self.process_group is not None:
            try:
                os.killpg(self.process_group, signum)
            except OSError:
                pass

    def _handle(self, signum, frame):
        del frame
        if self.received is None:
            self.received = signum
        self._forward(signum)

    def __enter__(self):
        for signum in _SIGNALS:
            self.previous[signum] = signal.signal(signum, self._handle)
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        del exc_type, exc_value, traceback
        if self._process_group_exists():
            self._terminate_process_group()
        for signum, handler in self.previous.items():
            signal.signal(signum, handler)

    def raise_if_cancelled(self):
        if self.received is not None:
            raise CheckM1Cancelled(self.received)

    def _terminate_process_group(self):
        if not self._process_group_exists():
            return
        signum = self.received or signal.SIGTERM
        self._forward(signum)
        deadline = time.monotonic() + CHILD_TERMINATION_GRACE_SECONDS
        while self._process_group_exists() and time.monotonic() < deadline:
            if self.process is not None:
                self.process.poll()
            time.sleep(0.05)
        if self._process_group_exists():
            try:
                os.killpg(self.process_group, signal.SIGKILL)
            except OSError:
                pass
        while self._process_group_exists():
            if self.process is not None:
                self.process.poll()
            time.sleep(0.05)

    def _popen_without_signal_gap(self, command, environ, stdout=None, stderr=None):
        old_mask = None
        restore_child_mask = None
        if hasattr(signal, "pthread_sigmask"):
            old_mask = signal.pthread_sigmask(signal.SIG_BLOCK, _SIGNALS)

            def restore_child_mask():
                signal.pthread_sigmask(signal.SIG_SETMASK, old_mask)

        try:
            self.raise_if_cancelled()
            process = subprocess.Popen(
                command,
                env=environ,
                start_new_session=True,
                stdout=stdout,
                stderr=stderr,
                preexec_fn=restore_child_mask,
            )
            self.process = process
            # start_new_session makes the child PID the stable process-group ID.
            # Keep it after the group leader exits so descendants remain protected.
            self.process_group = process.pid
        finally:
            if old_mask is not None:
                signal.pthread_sigmask(signal.SIG_SETMASK, old_mask)
        return process

    def _wait(self, process, capture):
        cancelled_at = None
        stdout = b""
        stderr = b""
        while True:
            try:
                if capture:
                    stdout, stderr = process.communicate(timeout=LOCK_POLL_SECONDS)
                    returncode = process.returncode
                else:
                    returncode = process.wait(timeout=LOCK_POLL_SECONDS)
                if not self._process_group_exists():
                    break
            except subprocess.TimeoutExpired:
                pass
            if self.received is not None:
                if cancelled_at is None:
                    cancelled_at = time.monotonic()
                    self._forward(self.received)
                elif time.monotonic() - cancelled_at >= CHILD_TERMINATION_GRACE_SECONDS:
                    try:
                        os.killpg(self.process_group, signal.SIGKILL)
                    except OSError:
                        pass
            if process.poll() is not None and self._process_group_exists():
                time.sleep(LOCK_POLL_SECONDS)
        if capture and self.received is not None:
            # communicate() after termination drains the small readiness probe pipes.
            stdout, stderr = process.communicate()
        return returncode, stdout, stderr

    def run(self, command, environ=None):
        if self.received is not None:
            return 128 + self.received
        process = self._popen_without_signal_gap(command, environ)
        returncode, _, _ = self._wait(process, False)
        self.process = None
        self.process_group = None
        if self.received is not None:
            return 128 + self.received
        return returncode

    def run_capture(self, command, environ=None):
        if self.received is not None:
            return 128 + self.received, "", ""
        process = self._popen_without_signal_gap(
            command, environ, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        returncode, stdout, stderr = self._wait(process, True)
        self.process = None
        self.process_group = None
        if self.received is not None:
            returncode = 128 + self.received
        return (
            returncode,
            stdout.decode("utf-8", "replace"),
            stderr.decode("utf-8", "replace"),
        )


def _native_database_probe(database, environ, forwarder):
    model = Path(database) / "hmms" / "checkm.hmm"
    try:
        first, _ = _bounded_edges(model)
    except OSError as error:
        raise CheckM1Error("Cannot inspect CheckM HMM for native probe: %s" % error)
    match = re.search(rb"^ACC\s+(\S+)", first, re.MULTILINE)
    if match is None:
        raise CheckM1Error("CheckM HMM has no accession for the bounded native probe.")
    key = match.group(1).decode("ascii", "strict")
    command = [_which("hmmfetch", environ), str(model), key]
    returncode, stdout, stderr = forwarder.run_capture(command, environ=environ)
    forwarder.raise_if_cancelled()
    if (
        returncode != 0
        or not stdout.startswith("HMMER3/")
        or not stdout.rstrip().endswith("//")
    ):
        raise CheckM1Error(
            "Native hmmfetch readiness probe failed for %s (exit %s): %s"
            % (database, returncode, stderr.strip())
        )


def _register_database(database, config_path, marker_root, forwarder, environ):
    """Atomically register, then verify fresh native defaults and marker parsing.

    CheckM 1.1.3's ``DBManager(set_path=...)`` rewrites the selected database's
    ``.dmanifest`` (twice in the pinned implementation).  Registration is the
    package-local ``DATA_CONFIG`` contract, so update that file atomically after
    readiness validation rather than mutating an already validated reference.
    """
    forwarder.raise_if_cancelled()
    try:
        values = json.loads(config_path.read_text(encoding="utf-8"))
        values["dataRoot"] = str(database)
        payload = json.dumps(values, sort_keys=True).encode("utf-8")
        _atomic_registration_write(config_path, payload, "register selected")
    except (OSError, ValueError, TypeError) as error:
        raise CheckM1Error(
            "CheckM DATA_CONFIG registration failed for %s: %s" % (database, error)
        )
    forwarder.raise_if_cancelled()

    try:
        configured = json.loads(config_path.read_text(encoding="utf-8"))["dataRoot"]
    except (OSError, ValueError, KeyError, TypeError) as error:
        raise CheckM1Error(
            "CheckM registration produced an unreadable DATA_CONFIG at %s: %s"
            % (config_path, error)
        )
    if _canonical(configured) != _canonical(database):
        raise CheckM1Error(
            "CheckM DATA_CONFIG did not retain the selected directory: expected %s, got %s"
            % (database, configured)
        )

    probe_code = (
        "from __future__ import print_function; "
        "import json,sys; "
        "from checkm.defaultValues import DefaultValues; "
        "from checkm.markerSets import MarkerSetParser; "
        "p=MarkerSetParser(); "
        "counts=[len(p.parseTaxonomicMarkerSetFile(x).getMarkerGenes()) for x in sys.argv[2:]]; "
        "print(json.dumps({'root':DefaultValues.CHECKM_DATA_DIR,'counts':counts}))"
    )
    probe_environment = dict(environ)
    probe_environment["PYTHONPATH"] = str(config_path.parent.parent)
    probe_environment["PYTHONNOUSERSITE"] = "1"
    probe = [
        sys.executable,
        "-c",
        probe_code,
        str(database),
        str(marker_root / "bacteria.ms"),
        str(marker_root / "archaea.ms"),
    ]
    returncode, stdout, stderr = forwarder.run_capture(probe, environ=probe_environment)
    forwarder.raise_if_cancelled()
    lines = [line.strip() for line in stdout.splitlines() if line.strip()]
    try:
        values = json.loads(lines[-1]) if lines else {}
        resolved = values["root"]
        counts = values["counts"]
    except (ValueError, KeyError, TypeError):
        resolved = ""
        counts = []
    if (
        returncode != 0
        or _canonical(resolved or ".") != _canonical(database)
        or len(counts) != 2
        or any(count < 40 for count in counts)
    ):
        raise CheckM1Error(
            "Fresh CheckM DefaultValues/marker verification failed for %s "
            "(exit %s, resolved %r, marker counts %r): %s"
            % (database, returncode, resolved, counts, stderr.strip())
        )


def _run_native_hook(
    hook,
    staged_payload,
    environment_prefix,
    config_path,
    snapshot,
    forwarder,
):
    """Run the installed legacy hook with package registration fully confined."""
    fake_prefix = Path(tempfile.mkdtemp(prefix="meew-checkm1-hook-", dir="/tmp"))
    alias = fake_prefix / "data"
    python_root = fake_prefix / "pythonpath"
    isolated_package = python_root / "checkm"
    alias.symlink_to(staged_payload, target_is_directory=True)
    python_root.mkdir()
    shutil.copytree(str(config_path.parent), str(isolated_package))
    (fake_prefix / "home").mkdir()
    hook_environment = os.environ.copy()
    hook_environment["PREFIX"] = str(fake_prefix)
    hook_environment["CHECKM_DATA_DIR"] = str(alias)
    hook_environment.pop("CHECKM_DATA_PATH", None)
    hook_environment["PYTHONPATH"] = os.pathsep.join(
        [str(python_root), hook_environment.get("PYTHONPATH", "")]
    ).rstrip(os.pathsep)
    hook_environment["PYTHONNOUSERSITE"] = "1"
    hook_environment["PYTHONDONTWRITEBYTECODE"] = "1"
    hook_environment["HOME"] = str(fake_prefix / "home")
    hook_environment["XDG_CONFIG_HOME"] = str(fake_prefix / "home" / ".config")
    hook_environment["PATH"] = os.pathsep.join(
        [str(Path(environment_prefix) / "bin"), hook_environment.get("PATH", "")]
    )
    remove_fake_prefix = True
    print("INFO: invoking installed CheckM Bioconda post-link setup hook: %s" % hook)
    print("INFO: native hook staging prefix: %s" % fake_prefix)
    try:
        returncode = forwarder.run(["bash", str(hook)], environ=hook_environment)
        messages = fake_prefix / ".messages.txt"
        if messages.is_file():
            print(
                "INFO: native CheckM hook diagnostics:\n%s"
                % messages.read_text(encoding="utf-8", errors="replace")
            )
        forwarder.raise_if_cancelled()
        if returncode != 0:
            raise CheckM1Error(
                "Installed CheckM post-link setup failed with exit status %s."
                % returncode
            )
        isolated_config = isolated_package / "DATA_CONFIG"
        try:
            isolated_root = json.loads(isolated_config.read_text(encoding="utf-8"))[
                "dataRoot"
            ]
        except (OSError, ValueError, KeyError, TypeError) as error:
            raise CheckM1Error(
                "Confined native hook registration is unreadable: %s" % error
            )
        if _canonical(isolated_root) != _canonical(staged_payload):
            raise CheckM1Error(
                "Confined native hook did not register its staging destination: %s"
                % isolated_root
            )
        if config_path.read_bytes() != snapshot:
            try:
                _restore_registration(config_path, snapshot)
            except BaseException as restore_error:
                remove_fake_prefix = False
                raise CheckM1RegistrationProtectionError(
                    "Installed hook unexpectedly changed the real environment "
                    "registration and rollback failed: %s. Confined prefix %s "
                    "was retained so DATA_CONFIG does not point to deleted staging."
                    % (restore_error, fake_prefix)
                )
            raise CheckM1Error(
                "Installed hook unexpectedly changed the real environment registration; "
                "the prior bytes were restored and setup was aborted."
            )
    finally:
        if remove_fake_prefix:
            shutil.rmtree(str(fake_prefix), ignore_errors=True)


def _destination_state(destination):
    if not destination.exists():
        return "missing"
    if not destination.is_dir():
        raise CheckM1Error(
            "CheckM 1 database destination is not a directory: %s" % destination
        )
    try:
        if not any(destination.iterdir()):
            return "empty"
    except OSError as error:
        raise CheckM1Error(
            "Cannot inspect CheckM 1 destination %s: %s" % (destination, error)
        )
    problems = database_problems(destination)
    if problems:
        raise CheckM1Error(
            "Refusing to merge, delete, or replace nonempty incomplete/corrupt "
            "CheckM 1 database at %s. Move it aside or select an empty exact "
            "database_path; readiness defect(s): %s"
            % (destination, "; ".join(problems))
        )
    return "valid"


def _prepare_database(
    destination,
    download_enabled,
    hook,
    environment_prefix,
    config_path,
    snapshot,
    forwarder,
    environ,
):
    forwarder.raise_if_cancelled()
    state = _destination_state(destination)
    if state == "valid":
        validate_database(destination)
        _native_database_probe(destination, environ, forwarder)
        print(
            "INFO: reusing complete CheckM 1 database without downloader: %s"
            % destination
        )
        return
    if not download_enabled:
        raise CheckM1Error(
            "CheckM 1 database is %s at %s and checkm1.download_database is false."
            % (state, destination)
        )

    forwarder.raise_if_cancelled()
    try:
        destination.parent.mkdir(parents=True, exist_ok=True)
        stage_root = Path(
            tempfile.mkdtemp(
                prefix=".%s.meew-stage-" % destination.name,
                dir=str(destination.parent),
            )
        )
    except OSError as error:
        raise CheckM1Error(
            "Cannot create protected CheckM 1 staging beside %s: %s"
            % (destination, error)
        )
    staged_payload = stage_root / "payload"
    staged_payload.mkdir()
    published = False
    preserve_staging = False
    try:
        _run_native_hook(
            hook,
            staged_payload,
            environment_prefix,
            config_path,
            snapshot,
            forwarder,
        )
        forwarder.raise_if_cancelled()
        validate_database(staged_payload)
        _native_database_probe(staged_payload, environ, forwarder)
        forwarder.raise_if_cancelled()

        try:
            os.replace(str(staged_payload), str(destination))
            published = True
        except OSError as error:
            raise CheckM1Error(
                "Validated CheckM 1 staging could not be atomically published to %s: %s"
                % (destination, error)
            )
        print(
            "INFO: atomically published validated CheckM 1 database: %s" % destination
        )
    except CheckM1RegistrationProtectionError:
        preserve_staging = True
        raise
    finally:
        if not published and not preserve_staging and staged_payload.exists():
            shutil.rmtree(str(staged_payload), ignore_errors=True)
        if not preserve_staging:
            shutil.rmtree(str(stage_root), ignore_errors=True)


def _validate_cleanup_path(
    cleanup_path, database, environment_prefix, protected_paths, tool_install_root
):
    if not cleanup_path:
        return None
    requested = Path(cleanup_path).expanduser().absolute()
    if requested.is_symlink():
        raise CheckM1Error(
            "MetaBinner cleanup path must not be a symlink: %s" % requested
        )
    cleanup = requested.resolve(strict=False)
    if cleanup == cleanup.parent:
        raise CheckM1Error(
            "MetaBinner cleanup path cannot be a filesystem root: %s" % cleanup
        )
    protected = [database, environment_prefix, tool_install_root]
    protected.extend(value for value in protected_paths if value)
    overlaps = [
        str(value) for value in protected if value and paths_overlap(cleanup, value)
    ]
    if overlaps:
        raise CheckM1Error(
            "Refusing MetaBinner output cleanup overlap at %s with protected path(s): %s"
            % (cleanup, ", ".join(overlaps))
        )
    if cleanup.exists() and not cleanup.is_dir():
        raise CheckM1Error("MetaBinner cleanup path is not a directory: %s" % cleanup)
    return cleanup


def prepare_and_run(
    database,
    download_enabled,
    protected_paths,
    command,
    environment_prefix=None,
    tool="comebin",
    tool_install_root="",
    native_paths=None,
    cleanup_path="",
):
    """Prepare/register CheckM and retain protection through one binner job."""
    prefix = _canonical(environment_prefix or sys.prefix)
    environment = os.environ.copy()
    with _SignalForwarder() as forwarder:
        forwarder.raise_if_cancelled()
        validate_native_paths(
            database,
            tool,
            tool_install_root=tool_install_root,
            native_paths=native_paths,
        )
        destination = validate_destination(database, protected_paths, prefix)
        cleanup = _validate_cleanup_path(
            cleanup_path,
            destination,
            prefix,
            protected_paths,
            tool_install_root,
        )
        installation = validate_native_installation(
            tool,
            command,
            environment,
            prefix,
            tool_install_root=tool_install_root,
        )
        forwarder.raise_if_cancelled()

        with _exclusive_lock(
            environment_lock_path(prefix), "environment", cancellation=forwarder
        ):
            forwarder.raise_if_cancelled()
            config_path = _checkm_data_config_path()
            _validate_registration_location(config_path, prefix)
            snapshot = _registration_snapshot(config_path)
            with _exclusive_lock(
                cache_lock_path(destination), "cache", cancellation=forwarder
            ):
                _prepare_database(
                    destination,
                    download_enabled,
                    installation["hook"],
                    prefix,
                    config_path,
                    snapshot,
                    forwarder,
                    environment,
                )

            forwarder.raise_if_cancelled()
            _isolate_registration_hardlinks(config_path, snapshot)
            try:
                _register_database(
                    destination,
                    config_path,
                    installation["marker_root"],
                    forwarder,
                    environment,
                )
            except BaseException:
                try:
                    _restore_registration(config_path, snapshot)
                except BaseException as restore_error:
                    raise CheckM1Error(
                        "CheckM setup failed and prior registration rollback also failed: %s"
                        % restore_error
                    )
                raise

            print(
                "INFO: CheckM 1 is ready and registered before %s starts: %s"
                % (tool, destination)
            )
            forwarder.raise_if_cancelled()
            if cleanup is not None and cleanup.exists():
                shutil.rmtree(str(cleanup))
            forwarder.raise_if_cancelled()
            return forwarder.run(command, environ=environment)


def _parse_boolean(value):
    normalized = value.strip().lower()
    if normalized in ("1", "true", "yes", "on"):
        return True
    if normalized in ("0", "false", "no", "off"):
        return False
    raise argparse.ArgumentTypeError("expected true or false")


def argument_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Safely prepare/register CheckM 1 and run one COMEBin or MetaBinner command."
        )
    )
    parser.add_argument("--database", required=True)
    parser.add_argument("--download-database", type=_parse_boolean, default=True)
    parser.add_argument("--tool", choices=SUPPORTED_TOOLS, default="comebin")
    parser.add_argument("--tool-install-root", default="")
    parser.add_argument("--protected-path", action="append", default=[])
    parser.add_argument("--native-path", action="append", default=[])
    parser.add_argument("--cleanup-path", default="")
    parser.add_argument("command", nargs=argparse.REMAINDER)
    return parser


def main(arguments=None):
    parsed = argument_parser().parse_args(arguments)
    command = list(parsed.command)
    if command and command[0] == "--":
        command.pop(0)
    try:
        return prepare_and_run(
            parsed.database,
            parsed.download_database,
            parsed.protected_path,
            command,
            tool=parsed.tool,
            tool_install_root=parsed.tool_install_root,
            native_paths=parsed.native_path,
            cleanup_path=parsed.cleanup_path,
        )
    except CheckM1Cancelled as error:
        print("ERROR: %s" % error, file=sys.stderr)
        return error.returncode
    except CheckM1Error as error:
        print("ERROR: %s" % error, file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
