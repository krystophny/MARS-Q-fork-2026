"""Behavioral provenance checks ported from accepted MARS-F commit 8c45256."""
import importlib.util
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
BUILD_SPEC = importlib.util.spec_from_file_location(
    "build_with_provenance", ROOT / "tools" / "build_with_provenance.py"
)
assert BUILD_SPEC is not None and BUILD_SPEC.loader is not None
build_provenance = importlib.util.module_from_spec(BUILD_SPEC)
BUILD_SPEC.loader.exec_module(build_provenance)


class BuildProvenance(unittest.TestCase):
    def test_build_manifest_binds_the_record_to_the_binary(self) -> None:
        with tempfile.TemporaryDirectory(prefix="mars-manifest-") as temporary:
            directory = Path(temporary)
            binary = directory / "marsq-test.x"
            binary.write_bytes(b"test MARS executable")
            manifest_path = directory / "marsq-test.x.provenance.json"
            manifest = {
                "schema": build_provenance.SCHEMA,
                "artifact": build_provenance.artifact_record(binary),
            }
            build_provenance.write_manifest(manifest_path, manifest)
            validated = build_provenance.validate_manifest(manifest_path)
            self.assertEqual(
                validated["artifact"]["sha256"],
                build_provenance.sha256_file(binary),
            )
            binary.write_bytes(b"fake MARS executable")
            with self.assertRaisesRegex(ValueError, "mismatch"):
                build_provenance.validate_manifest(manifest_path)

    def test_linked_libraries_use_otool_on_darwin(self) -> None:
        binary = Path("/tmp/marsq-test.x")
        with (
            mock.patch.object(
                build_provenance.platform, "system", return_value="Darwin"
            ),
            mock.patch.object(
                build_provenance.shutil, "which", return_value="/usr/bin/otool"
            ),
            mock.patch.object(
                build_provenance,
                "run_output",
                return_value=(
                    "/tmp/marsq-test.x:\n"
                    "\t/usr/lib/libSystem.B.dylib (compatibility version 1.0.0)"
                ),
            ) as run_output,
        ):
            command, libraries = build_provenance.linked_libraries(binary)

        self.assertEqual(command, ["/usr/bin/otool", "-L", str(binary)])
        self.assertEqual(len(libraries), 2)
        run_output.assert_called_once_with(command, required=False)

    def test_linked_libraries_use_ldd_on_linux(self) -> None:
        binary = Path("/tmp/marsq-test.x")
        with (
            mock.patch.object(
                build_provenance.platform, "system", return_value="Linux"
            ),
            mock.patch.object(
                build_provenance.shutil, "which", return_value="/usr/bin/ldd"
            ),
            mock.patch.object(
                build_provenance, "run_output", return_value="libgfortran.so.5"
            ) as run_output,
        ):
            command, libraries = build_provenance.linked_libraries(binary)

        self.assertEqual(command, ["/usr/bin/ldd", str(binary)])
        self.assertEqual(libraries, ["libgfortran.so.5"])
        run_output.assert_called_once_with(command, required=False)

    def test_linked_libraries_allow_missing_inspection_tool(self) -> None:
        with mock.patch.object(build_provenance.shutil, "which", return_value=None):
            command, libraries = build_provenance.linked_libraries(Path("marsq.x"))

        self.assertIsNone(command)
        self.assertEqual(libraries, [])

    def test_mpi_wrapper_path_preserves_argv_zero_dispatch(self) -> None:
        wrapper = shutil.which("mpif90")
        if wrapper is None:
            self.skipTest("mpif90 is not available")
        self.assertEqual(
            build_provenance.resolve_program("mpif90"), Path(wrapper).absolute()
        )
