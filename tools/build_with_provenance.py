#!/usr/bin/env python3
"""Build a MARS executable and bind it to machine-readable provenance."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
SOURCE_DIR = ROOT / "MarsQ_2FK"
SCHEMA = "mars-build-provenance-v1"
PROFILES = {
    "gnu": {
        "binary": "marsq-gnu.x",
        "compiler": "mpif90",
        "flags": [
            "-O2",
            "-fdefault-real-8",
            "-fdefault-double-8",
            "-fopenmp",
            "-ffixed-line-length-none",
            "-ffree-line-length-none",
            "-fallow-argument-mismatch",
        ],
    },
    "gnu-debug": {
        "binary": "marsq-gnu-debug.x",
        "compiler": "mpif90",
        "flags": [
            "-O0",
            "-g",
            "-fdefault-real-8",
            "-fdefault-double-8",
            "-fopenmp",
            "-ffixed-line-length-none",
            "-ffree-line-length-none",
            "-fallow-argument-mismatch",
            "-fcheck=all",
            "-fbacktrace",
        ],
    },
    "ifx": {
        "binary": "marsq-ifx.x",
        "compiler": "mpiifx",
        "flags": ["-O2", "-real-size", "64", "-qopenmp", "-extend-source", "132"],
    },
    "nvhpc": {
        "binary": "marsq-nvhpc.x",
        "compiler": "mpifort",
        "flags": ["-O1", "-r8", "-mp", "-Mextend"],
    },
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run_output(
    command: list[str], *, cwd: Path | None = None, required: bool = True
) -> str | None:
    result = subprocess.run(
        command,
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if result.returncode != 0:
        if required:
            rendered = shlex.join(command)
            raise RuntimeError(
                f"command failed ({result.returncode}): {rendered}\n{result.stdout}"
            )
        return None
    return result.stdout.strip()


def resolve_program(name: str) -> Path:
    resolved = shutil.which(name)
    if resolved is None:
        raise RuntimeError(f"required build program is not on PATH: {name}")
    # MPI compiler wrappers select their configuration from argv[0].  Keep the
    # named wrapper path here; invoking its resolved opal_wrapper symlink breaks
    # that dispatch.  The real path is recorded separately in the manifest.
    return Path(resolved).absolute()


def git_state() -> dict[str, Any]:
    status_text = run_output(
        ["git", "status", "--porcelain=v1", "--untracked-files=all"], cwd=ROOT
    )
    assert status_text is not None
    diff = subprocess.run(
        ["git", "diff", "--binary", "HEAD", "--"],
        cwd=ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    ).stdout
    untracked_text = run_output(
        ["git", "ls-files", "--others", "--exclude-standard"], cwd=ROOT
    )
    assert untracked_text is not None
    untracked = []
    for relative in untracked_text.splitlines():
        path = ROOT / relative
        entry: dict[str, Any] = {"path": relative}
        if path.is_symlink():
            entry["symlink_target"] = os.readlink(path)
        elif path.is_file():
            entry["bytes"] = path.stat().st_size
            entry["sha256"] = sha256_file(path)
        untracked.append(entry)
    remotes = {}
    remote_text = run_output(["git", "remote"], cwd=ROOT)
    assert remote_text is not None
    for remote in remote_text.splitlines():
        url = run_output(["git", "remote", "get-url", remote], cwd=ROOT)
        assert url is not None
        remotes[remote] = url
    commit = run_output(["git", "rev-parse", "HEAD"], cwd=ROOT)
    branch = run_output(
        ["git", "symbolic-ref", "--quiet", "--short", "HEAD"],
        cwd=ROOT,
        required=False,
    )
    describe = run_output(
        ["git", "describe", "--always", "--dirty", "--tags"], cwd=ROOT
    )
    assert commit is not None and describe is not None
    return {
        "repository": str(ROOT),
        "commit": commit,
        "branch": branch,
        "describe": describe,
        "dirty": bool(status_text),
        "status_porcelain": status_text.splitlines(),
        "diff_head_sha256": hashlib.sha256(diff).hexdigest(),
        "untracked": untracked,
        "remotes": remotes,
    }


def source_fingerprint(state: dict[str, Any]) -> tuple[Any, ...]:
    return (
        state["commit"],
        state["dirty"],
        tuple(state["status_porcelain"]),
        state["diff_head_sha256"],
        json.dumps(state["untracked"], sort_keys=True),
    )


def compiler_record(compiler: Path) -> dict[str, Any]:
    version = run_output([str(compiler), "--version"])
    assert version is not None
    record: dict[str, Any] = {
        "path": str(compiler),
        "realpath": str(compiler.resolve()),
        "version": version,
    }
    for option, key in (
        ("--showme:version", "wrapper_version"),
        ("--showme:command", "backend_command"),
        ("--showme:compile", "wrapper_compile_flags"),
        ("--showme:link", "wrapper_link_flags"),
    ):
        value = run_output([str(compiler), option], required=False)
        if value is not None:
            record[key] = value
    backend = record.get("backend_command")
    if backend:
        backend_words = shlex.split(backend)
        if backend_words:
            backend_path = shutil.which(backend_words[0])
            if backend_path is not None:
                record["backend_path"] = str(Path(backend_path).absolute())
                record["backend_realpath"] = str(Path(backend_path).resolve())
                backend_version = run_output(
                    [str(Path(backend_path).absolute()), "--version"], required=False
                )
                if backend_version is not None:
                    record["backend_version"] = backend_version
    return record


def artifact_record(binary: Path) -> dict[str, Any]:
    file_program = shutil.which("file")
    description = None
    if file_program is not None:
        description = run_output([file_program, "-b", str(binary)])
    return {
        "path": str(binary.resolve()),
        "bytes": binary.stat().st_size,
        "sha256": sha256_file(binary),
        "file_description": description,
    }


def linked_libraries(binary: Path) -> tuple[list[str] | None, list[str]]:
    if platform.system() == "Darwin":
        program_name = "otool"
        arguments = ["-L", str(binary)]
    else:
        program_name = "ldd"
        arguments = [str(binary)]

    program = shutil.which(program_name)
    if program is None:
        return None, []

    command = [program, *arguments]
    output = run_output(command, required=False)
    return command, output.splitlines() if output else []


def write_manifest(path: Path, manifest: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
        delete=False,
    ) as stream:
        temporary = Path(stream.name)
        json.dump(manifest, stream, indent=2, sort_keys=True)
        stream.write("\n")
    os.replace(temporary, path)


def validate_manifest(path: Path) -> dict[str, Any]:
    manifest = json.loads(path.read_text(encoding="utf-8"))
    if manifest.get("schema") != SCHEMA:
        raise ValueError(f"unsupported provenance schema in {path}")
    artifact = manifest.get("artifact")
    if not isinstance(artifact, dict):
        raise ValueError(f"missing artifact record in {path}")
    binary = Path(artifact.get("path", ""))
    if not binary.is_file():
        raise ValueError(f"provenance artifact does not exist: {binary}")
    actual_bytes = binary.stat().st_size
    if actual_bytes != artifact.get("bytes"):
        raise ValueError(
            f"artifact size mismatch: manifest={artifact.get('bytes')} actual={actual_bytes}"
        )
    actual_hash = sha256_file(binary)
    if actual_hash != artifact.get("sha256"):
        raise ValueError(
            f"artifact SHA-256 mismatch: manifest={artifact.get('sha256')} "
            f"actual={actual_hash}"
        )
    return manifest


def build(profile_name: str, build_dir: Path, require_clean: bool, *,
          pams_csc: bool = False, fortsparse_prefix: Path | None = None,
          fortsparse_libs: str = "") -> Path:
    if pams_csc and (profile_name not in ("gnu", "gnu-debug") or fortsparse_prefix is None):
        raise ValueError("PAMS CSC requires a GNU profile and --fortsparse-prefix")
    if not pams_csc and (fortsparse_prefix is not None or fortsparse_libs):
        raise ValueError("fortsparse options require --pams-csc")
    dependency = {"enabled": pams_csc}
    if pams_csc:
        assert fortsparse_prefix is not None
        fortsparse_prefix = fortsparse_prefix.expanduser().resolve()
        dependency.update(prefix=str(fortsparse_prefix), extra_link_flags=shlex.split(fortsparse_libs))
        dependency_files = [fortsparse_prefix / "lib/libfortsparse.a",
                            fortsparse_prefix / "bin/fortsparse_umfpack_helper"]
        dependency_files += sorted((fortsparse_prefix / "include").glob("*.mod"))
        if not dependency_files[2:]:
            raise ValueError("fortsparse prefix has no Fortran modules")
        dependency["files_sha256"] = {str(p): sha256_file(p) for p in dependency_files}
    profile = PROFILES[profile_name]
    compiler = resolve_program(profile["compiler"])
    make = resolve_program("make")
    python = Path(sys.executable).resolve()
    build_dir.mkdir(parents=True, exist_ok=True)
    binary = (build_dir / profile["binary"]).resolve()
    manifest_path = binary.with_name(f"{binary.name}.provenance.json")
    temporary_binary = build_dir / f".{binary.name}.{os.getpid()}.tmp"
    temporary_binary.unlink(missing_ok=True)

    source_before = git_state()
    if require_clean and source_before["dirty"]:
        raise RuntimeError("refusing production build from a dirty Git worktree")

    flags = list(profile["flags"])
    flags_string = " ".join(flags)
    clean_command = [str(make), "-C", str(SOURCE_DIR), "clean-objects"]
    build_command = [
        str(make),
        "-C",
        str(SOURCE_DIR),
        "--no-print-directory",
        f"PROG={temporary_binary}",
        f"F95={compiler}",
        f"F95FLAGS={flags_string}",
    ]
    build_command.append("PAMS_CSC=1" if pams_csc else "PAMS_CSC=0")
    if pams_csc:
        build_command += [f"FORTSPARSE_PREFIX={fortsparse_prefix}",
                          f"FORTSPARSE_LIBS={fortsparse_libs}"]
    wrapper_command = [
        str(python),
        str(Path(__file__).resolve()),
        "--profile",
        profile_name,
        "--build-dir",
        str(build_dir),
    ]
    if pams_csc:
        wrapper_command += ["--pams-csc", "--fortsparse-prefix", str(fortsparse_prefix),
                            "--fortsparse-libs=" + fortsparse_libs]
    if require_clean:
        wrapper_command.append("--require-clean")

    started_at = utc_now()
    start = time.monotonic()
    try:
        subprocess.run(clean_command, check=True)
        subprocess.run(build_command, check=True)
        if not temporary_binary.is_file():
            raise RuntimeError(
                f"build did not create expected executable: {temporary_binary}"
            )
        if pams_csc:
            if dependency["files_sha256"] != {str(p): sha256_file(p) for p in dependency_files}:
                raise RuntimeError("fortsparse modules/library/helper changed during build")
        source_after = git_state()
        if source_fingerprint(source_before) != source_fingerprint(source_after):
            raise RuntimeError("Git source state changed while MARS was compiling")
        os.replace(temporary_binary, binary)
    finally:
        temporary_binary.unlink(missing_ok=True)

    make_version = run_output([str(make), "--version"])
    assert make_version is not None
    library_command, library_lines = linked_libraries(binary)
    environment_names = (
        "PATH",
        "LD_LIBRARY_PATH",
        "LIBRARY_PATH",
        "CPATH",
        "FC",
        "F90",
        "OMPI_FC",
        "OMPI_FCFLAGS",
        "I_MPI_FC",
        "I_MPI_F90",
        "NVHPC",
        "ONEAPI_ROOT",
    )
    manifest = {
        "schema": SCHEMA,
        "artifact": artifact_record(binary),
        "build": {
            "profile": profile_name,
            "pams_csc": dependency,
            "flags": flags,
            "flags_string": flags_string,
            "compiler": compiler_record(compiler),
            "make": {"path": str(make), "version": make_version},
            "executed_commands": [clean_command, build_command],
            "wrapper_command": wrapper_command,
            "linked_libraries_command": library_command,
            "linked_libraries": library_lines,
            "environment": {
                name: os.environ[name]
                for name in environment_names
                if name in os.environ
            },
            "started_at_utc": started_at,
            "completed_at_utc": utc_now(),
            "elapsed_seconds": time.monotonic() - start,
        },
        "source": source_before,
        "source_unchanged_during_build": True,
        "host": {
            "hostname": platform.node(),
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python": {"path": str(python), "version": platform.python_version()},
        },
    }
    write_manifest(manifest_path, manifest)
    validate_manifest(manifest_path)
    print(f"built {binary}")
    print(f"provenance {manifest_path}")
    print(f"sha256 {manifest['artifact']['sha256']}")
    return manifest_path


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--profile", choices=sorted(PROFILES))
    action.add_argument(
        "--verify-manifest",
        type=Path,
        help="verify that a manifest still matches its recorded executable",
    )
    parser.add_argument(
        "--build-dir",
        type=Path,
        default=ROOT / "build",
        help="output directory (default: repository build directory)",
    )
    parser.add_argument(
        "--require-clean",
        action="store_true",
        help="refuse to build when tracked or untracked Git changes are present",
    )
    parser.add_argument("--pams-csc", action="store_true",
                        help="enable globally pivoted steady-response backend")
    parser.add_argument("--fortsparse-prefix", type=Path)
    parser.add_argument("--fortsparse-libs", default="",
                        help="additional native static-link flags")
    return parser.parse_args()


def main() -> int:
    arguments = parse_arguments()
    if arguments.verify_manifest is not None:
        manifest = validate_manifest(arguments.verify_manifest.resolve())
        print(
            f"verified {manifest['artifact']['path']} "
            f"sha256={manifest['artifact']['sha256']}"
        )
        return 0
    build_dir = arguments.build_dir.expanduser()
    if not build_dir.is_absolute():
        build_dir = (Path.cwd() / build_dir).resolve()
    assert arguments.profile is not None
    build(arguments.profile, build_dir, arguments.require_clean,
          pams_csc=arguments.pams_csc, fortsparse_prefix=arguments.fortsparse_prefix,
          fortsparse_libs=arguments.fortsparse_libs)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, ValueError, subprocess.SubprocessError) as error:
        print(f"build provenance error: {error}", file=sys.stderr)
        raise SystemExit(2) from error
