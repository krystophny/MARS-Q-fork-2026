"""Execute native KDIAMAG against physical thermal-rate and legacy-grid oracles."""
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


def production_routine(source, name):
    start = re.search(r"^      SUBROUTINE " + name + r"(?:\(|\s*$)", source, re.M).start()
    end = re.search(r"^      END(?: SUBROUTINE " + name + r")?\s*$", source[start:], re.M)
    return source[start:start + end.end()] + "\n"


class ThermalSpeciesCollision(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = Path(os.environ.get("MARS_ORACLE_SOURCE", ROOT)) / "MarsQ_2FK"
        scratch = ROOT / "build"
        scratch.mkdir(exist_ok=True)
        cls.temporary = tempfile.TemporaryDirectory(dir=scratch, prefix="thermal-rate-")
        cls.addClassCleanup(cls.temporary.cleanup)
        cls.work = Path(cls.temporary.name)
        source = (cls.source / "kinetic.f").read_text()
        kernel = production_routine(source, "KDIAMAG")
        if "SUBROUTINE KVALIDATESPECIESNU" in source:
            kernel += production_routine(source, "KVALIDATESPECIESNU")
        # On old sources the input metadata are local unused fixture variables;
        # the executed production KDIAMAG remains byte-for-byte unchanged.
        old_metadata = ""
        if "KNUSPECIES" not in (cls.source / "globalm.f").read_text():
            old_metadata = "integer :: KNUSPECIES\nreal*8 :: NUEFFSA(10)\n"
        default_profiles = ""
        if "KPROFSPEC" in (cls.source / "globalm.f").read_text():
            default_profiles = "KPROFSPEC=0"
        readback = ""
        if not old_metadata:
            kernel += production_routine((cls.source / "marsq.f").read_text(), "CALC_NUSTAR")
            readback = "if (KNUSPECIES==1) call CALC_NUSTAR"
        (cls.work / "kernel.f").write_text(kernel)
        (cls.work / "driver.f90").write_text(
            DRIVER.replace("!INPUT_METADATA", old_metadata).replace("!READBACK", readback)
                  .replace("!DEFAULT_PROFILES", default_profiles))
        (cls.work / "mpi_fixture.f90").write_text(
            "module MPIENV\ninteger :: ISMPIRUN,RANK\ninteger,parameter :: ROOT=0\nend module\n")
        modules = ["globalm.f", "dimensim.f", "rcomdm.f", "kineticm.f", "anisotropicm.f", "torquem.f"]
        command = [os.environ.get("FC", "gfortran"), "-O0", "-fcheck=all",
                   "-fdefault-real-8", "-fdefault-double-8", "-ffixed-line-length-none",
                   "-ffree-line-length-none", "-I" + str(cls.source)]
        command += [str(cls.source / name) for name in modules]
        command += ["mpi_fixture.f90", "kernel.f", "driver.f90", "-o", "oracle"]
        result = subprocess.run(command, cwd=cls.work, capture_output=True, text=True, timeout=60)
        if result.returncode:
            raise RuntimeError(result.stdout + result.stderr)

    def execute(self, opt_in, species=3, model=0, rates=(1, 2, 3), inutype=0, prof=0, nexp=0, kinetic=1, kntv=21):
        inputs = [opt_in, species, model, inutype, prof, nexp, kinetic, kntv, *rates]
        for name in ("PROFNUSTAR.OUT", "THERMAL_COLLISION_RATES.OUT"):
            (self.work / name).unlink(missing_ok=True)
        result = subprocess.run([str(self.work / "oracle")], cwd=self.work,
                                input=" ".join(map(str, inputs)) + "\n",
                                capture_output=True, text=True, timeout=5)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return result.stdout, result.stderr

    def test_native_physical_widths_and_legacy_profiles(self):
        # Independent NRL self-collision coefficient, with T expressed in eV.
        # This supplies effective rates, without asserting all-pair FP physics.
        charge = 1.602176634e-19
        epsilon0 = 8.8541878128e-12
        proton_mass = 1.67262192369e-27
        prefactor = 4 * math.sqrt(math.pi) / 3 * charge**2.5 / (4 * math.pi * epsilon0)**2
        def rate(mass, density=1e19, temperature=1000, z=1):
            return prefactor * density * 15 * z**4 / math.sqrt(proton_mass * mass) / temperature**1.5
        rates = (rate(2), rate(1 / 1836.15267343, 2e19, 1500), rate(3))
        self.assertAlmostEqual(rates[2] / rates[0], math.sqrt(2 / 3), places=14)
        for model in (0, 1):
            for opt_in, species, supplied in ((0, 2, rates), (1, 3, rates), (1, 3, (0, 0, 0))):
                with self.subTest(model=model, opt_in=opt_in, species=species, rates=supplied):
                    out, err = self.execute(opt_in, species, model, supplied)
                    self.assertNotIn("STOP", err)
                    rows = [list(map(float, line.split())) for line in out.splitlines()]
                    self.assertEqual(len(rows), 6)
                    for row in rows:
                        j, grid = map(int, row[:2])
                        # Integer axis is copied from the first half-grid value.
                        effective_grid = 2 if j == 1 and grid == 1 else grid
                        coordinate = ([.2, .6, 1] if effective_grid == 1 else [.4, .8, 1.2])[j - 1]
                        factor = 1 if model == 0 else 3 / (2 * coordinate)
                        if opt_in:
                            expected = [value * 5e-6 * factor for value in supplied]
                        else:
                            ion = ([.01, .02, .03] if effective_grid == 1 else [.015, .025, .035])[j - 1]
                            electron = ([.1, .2, .3] if effective_grid == 1 else [.15, .25, .35])[j - 1]
                            expected = [ion * factor, electron * factor, 0]
                        for index in (2, 0, 1):
                            actual, wanted = row[2 + index], expected[index]
                            self.assertAlmostEqual(actual, wanted, delta=max(1e-15, abs(wanted) * 2e-15))
                    if opt_in:
                        rates_file = self.work / "THERMAL_COLLISION_RATES.OUT"
                        self.assertTrue(rates_file.is_file(), "active rate readback missing")
                        actual_rates = [list(map(float, line.split())) for line in
                                        rates_file.read_text().splitlines() if not line.startswith("%")]
                        self.assertEqual(len(actual_rates), 3)
                        for index, row in enumerate(actual_rates):
                            self.assertEqual(int(row[0]), index + 1)
                            self.assertAlmostEqual(row[1], supplied[index], delta=max(1e-15, abs(supplied[index]) * 2e-15))
                            self.assertAlmostEqual(row[2], supplied[index] * 5e-6, delta=max(1e-15, abs(supplied[index] * 5e-6) * 2e-15))
                        diagnostic = (self.work / "PROFNUSTAR.OUT").read_text()
                        self.assertIn("DIAGNOSTIC OMITTED", diagnostic)
                        self.assertTrue(all(line.startswith("%") for line in diagnostic.splitlines()))

    def test_invalid_rate_and_profile_inputs_reject(self):
        cases = [(dict(rates=(-1, 2, 3)), "NONNEGATIVE"),
                 (dict(rates=("NaN", 2, 3)), "FINITE"),
                 (dict(rates=("Infinity", 2, 3)), "FINITE"),
                 (dict(rates=(1, 2, 0), inutype=1), "POSITIVE THERMAL"),
                 (dict(prof=4), "NPROFIE=0 AND NEXPV=0"),
                 (dict(nexp=1), "NPROFIE=0 AND NEXPV=0"),
                 (dict(kntv=11), "INCKIN=1 AND KNTV=21"),
                 (dict(kinetic=0), "INCKIN=1 AND KNTV=21"),
                 (dict(rates=(1e-320, 2, 3), inutype=1), "POSITIVE ACTIVE NUEFF")]
        for arguments, diagnostic in cases:
            with self.subTest(arguments=arguments):
                out, err = self.execute(1, **arguments)
                self.assertEqual(out, "")
                self.assertIn(diagnostic, err)


DRIVER = """
program thermal_rate_oracle
 use rcomdm
 use dimensim
 use globalm
 use kineticm
 use mpienv
 use anisotropicm
 implicit none
 integer :: j,k
!INPUT_METADATA
 NR=2;NRP1=3;ISMPIRUN=1;RANK=ROOT
!DEFAULT_PROFILES
 read(*,*) KNUSPECIES,NSPECIES,NUMODEL,INUTYPE,NPROFIE,NEXPV,INCKIN,KNTV,(NUEFFSA(k),k=1,3)
 allocate(CS(3),CSM(3),GNUI(3),GNUIM(3),GNUE(3),GNUEM(3))
 allocate(NUEFF(3,2,3),FFF(3,2),DFFF(3,2))
 allocate(EPSLONCA(3,2,3),EPSALPHA(3,2,3),AAK(3,2,3))
 allocate(DEPSALPHADPSI(3,2,3),DEPSCDPSI2(3,2,3),DEPSALPHADPSI2(3,2,3))
 allocate(OMEGASN(3,2,3),OMEGAST(3,2,3),OMEGASNA(3,2,3),OMEGASAA(3,2,3),OMEGASCA(3,2,3))
 allocate(ESPECIES_DEN(3,2,3),ESPECIES_TEM(3,2,3))
 CS=[.2d0,.6d0,1d0];CSM=[.4d0,.8d0,1.2d0];ASPCT=3d0;ZTAUA0=5d-6
 GNUI=[.01d0,.02d0,.03d0];GNUIM=[.015d0,.025d0,.035d0]
 GNUE=[.1d0,.2d0,.3d0];GNUEM=[.15d0,.25d0,.35d0]
 EPSLONCA=1d0;EPSALPHA=1d0;AAK=1d0
 ESPECIES_DEN=1d0;ESPECIES_DEN(:,:,2)=2d0;ESPECIES_TEM=1d0
 ESPECIES_M=2d0;ESPECIES_M(2)=1d0/1836.15267343d0;ESPECIES_M(3)=3d0
 ESPECIES_Z=1d0;ESPECIES_Z(2)=-1d0
 ISPECIES_F0=0;B0K=1d0;OMEGACI0=50d0;B0EXP=1d0;R0EXP=1d0
 call KDIAMAG
!READBACK
 do j=1,3
  do k=1,2
   print '(2I3,3ES24.16)',j,k,NUEFF(j,k,:)
  enddo
 enddo
end program
subroutine DFFFDPSI(mode)
 use globalm,only:DFFF
 integer :: mode
 DFFF=0d0
end subroutine
subroutine KWRITEEXPLICITPROFILES
 error stop 'unexpected explicit profile observer in default-off width test'
end subroutine
"""

if __name__ == "__main__":
    unittest.main()
