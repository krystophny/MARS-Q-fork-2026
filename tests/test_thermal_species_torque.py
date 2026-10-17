"""Execute native KNTV=21 aggregation against manufactured species work.

The unchanged CALCDWKCOMP block selects Hermitian work at zero orbit width
and dissipative-kernel work for either finite-orbit-width flag.  Three
independent species responses expose truncation after the historical ion and
electron slots.  Only Maxwellian species contribute to the thermal total;
the original first-two-slot outputs retain their separate responses.
"""
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
HERMITIAN = (1.0, 2.0, 3.0)
DISSIPATIVE = (4.0, 5.0, 6.0)
WEAK = (7.0, 8.0, 9.0)
SURFACES = (0.25, 0.75)
CELL_SCALE = (1.0, 10.0)


def torque_source(source):
    """Copy native array declarations and aggregation, including its closing IF."""
    text = (source / "MarsQ_2FK" / "kinetic.f").read_text()
    routine = re.search(r"^      SUBROUTINE CALCDWKCOMP\(", text, re.M)
    if routine is None:
        raise AssertionError("native CALCDWKCOMP routine not found")
    end = text.index("      END SUBROUTINE CALCDWKCOMP", routine.end())
    body = text[routine.start():end]
    declaration = re.search(r"^      REAL\*8\s+TH\([^\n]+$", body, re.M)
    if declaration is None:
        raise AssertionError("native torque work-array declaration not found")
    start = re.search(r"^      TNORM\s*=\s*-1\.\*RNTOR/\(4\.\*PI\*PI\)\s*$",
                      body, re.M)
    if start is None:
        raise AssertionError("native torque normalization not found")
    section = re.search(r"^C\s+(?:SAME-DRIVE WORK|OUTPUT PROFILES)\b",
                        body[start.end():], re.M)
    if section is None:
        raise AssertionError("native torque section terminator not found")
    return declaration.group(), body[start.start():start.end() + section.start()]


DRIVER = r"""
      PROGRAM THERMAL_SPECIES_TORQUE
      IMPLICIT NONE
      INTEGER, PARAMETER :: NSPECIES=3,NR=2
      REAL*8, PARAMETER :: PI=ACOS(-1.0D0)
      INTEGER KP,IS,I,INDX,IH,NSP,FID,FIDS,IFOWP,IFOWT
      INTEGER ISPECIES_F0(NSPECIES),INDXDWKC(NSPECIES,10)
      INTEGER INDXHWKC(NSPECIES,5),ASSIGNFREEFILEUNIT
@DECLARATION@
      REAL*8 TORQUENTVI(NR),TORQUENTVE(NR),TORQUENTV(NR)
      REAL*8 CSM(NR),RNTOR,FACTOR,HW,DW,WW,POLARITY
      COMPLEX*16 DWPPARY(NR,30),DWPPERY(NR,30),DWHERM(NR,3)
      CHARACTER*32 ARG
      CALL GET_COMMAND_ARGUMENT(1,ARG)
      READ(ARG,*) IFOWP
      CALL GET_COMMAND_ARGUMENT(2,ARG)
      READ(ARG,*) IFOWT
      CALL GET_COMMAND_ARGUMENT(3,ARG)
      ISPECIES_F0=0
      READ(ARG,*) ISPECIES_F0(3)
      CALL GET_COMMAND_ARGUMENT(4,ARG)
      READ(ARG,*) POLARITY
C     TNORM=POLARITY EXACTLY: THIS FIXTURE TESTS AGGREGATION, NOT THE
C     TORQUE NORMALIZATION ALREADY COVERED BY THE REAL-FIELD ORACLE.
      RNTOR=-4.0D0*PI*PI*POLARITY
      CSM=(/0.25D0,0.75D0/)
      INDXDWKC=-1
      INDXHWKC=-1
      DWPPARY=0
      DWPPERY=0
      DWHERM=0
      TORQUENTVI=-999
      TORQUENTVE=-999
      TORQUENTV=-999
      DO KP=1,NSPECIES
         INDXDWKC(KP,3)=10*(KP-1)+3
         INDXDWKC(KP,8)=10*(KP-1)+8
         INDXHWKC(KP,3)=KP
         DO IS=1,NR
            FACTOR=10.0D0**(IS-1)
            HW=FACTOR*DBLE(KP)
            DW=FACTOR*DBLE(KP+3)
            WW=FACTOR*DBLE(KP+6)
            DWHERM(IS,KP)=CMPLX(11.0D0+KP,-HW,KIND=8)
            INDX=INDXDWKC(KP,3)
            DWPPARY(IS,INDX)=CMPLX(100.0D0+KP,-0.25D0*WW,KIND=8)
            DWPPERY(IS,INDX)=CMPLX(-35.0D0,-0.75D0*WW,KIND=8)
            INDX=INDXDWKC(KP,8)
            DWPPARY(IS,INDX)=CMPLX(200.0D0+KP,-0.25D0*DW,KIND=8)
            DWPPERY(IS,INDX)=CMPLX(-70.0D0,-0.75D0*DW,KIND=8)
         ENDDO
      ENDDO
C     THE NATIVE BLOCK STARTS INSIDE THE KNTV=21 CONDITIONAL.
      IF (.TRUE.) THEN
@PRODUCTION@
      DO IS=1,NR
         WRITE(*,'(A,I4,3ES26.17E3)') 'TORQUE ',IS,
     &      TORQUENTVI(IS),TORQUENTVE(IS),TORQUENTV(IS)
      ENDDO
      END PROGRAM THERMAL_SPECIES_TORQUE

      INTEGER FUNCTION ASSIGNFREEFILEUNIT()
      IMPLICIT NONE
      INTEGER U
      LOGICAL OPENED
      DO U=20,99
         INQUIRE(UNIT=U,OPENED=OPENED)
         IF (.NOT.OPENED) THEN
            ASSIGNFREEFILEUNIT=U
            RETURN
         ENDIF
      ENDDO
      STOP 'NO FREE FIXTURE UNIT'
      END FUNCTION ASSIGNFREEFILEUNIT
"""


def rows(path):
    if not path.exists():
        return None
    return [[float(x) for x in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("%")]


@unittest.skipUnless(shutil.which("gfortran"), "native fixture requires gfortran")
class ThermalSpeciesTorque(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        source = Path(os.environ.get("MARS_ORACLE_SOURCE", ROOT)).resolve()
        if source.name == "MarsQ_2FK":
            source = source.parent
        scratch = Path(os.environ.get("TMPDIR", ROOT / "build")).resolve()
        if scratch == Path("/tmp") or Path("/tmp") in scratch.parents:
            raise AssertionError("native fixture scratch must be outside /tmp")
        scratch.mkdir(parents=True, exist_ok=True)
        cls.temporary = tempfile.TemporaryDirectory(prefix="thermal-torque-", dir=scratch)
        cls.addClassCleanup(cls.temporary.cleanup)
        cls.work = Path(cls.temporary.name)
        declaration, block = torque_source(source)
        (cls.work / "oracle.f").write_text(
            DRIVER.replace("@DECLARATION@", declaration).replace("@PRODUCTION@", block))
        command = ["gfortran", "-O0", "-fcheck=all", "-std=legacy",
                   "-fdefault-real-8", "-fdefault-double-8",
                   "-ffixed-line-length-none", "-ffree-line-length-none",
                   "-fallow-argument-mismatch", "oracle.f", "-o", "oracle"]
        result = subprocess.run(command, cwd=cls.work, capture_output=True,
                                text=True, timeout=120)
        if result.returncode:
            raise AssertionError(result.stdout + result.stderr)

    def run_case(self, ifowp, ifowt, third_f0=0, polarity=1):
        for name in ("NTV_TORQUE_FORMS.OUT", "NTV_TORQUE_SPECIES.OUT"):
            (self.work / name).unlink(missing_ok=True)
        result = subprocess.run([str(self.work / "oracle"), str(ifowp), str(ifowt),
                                 str(third_f0), str(polarity)], cwd=self.work,
                                capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        native = [[float(x) for x in line.split()[1:]]
                  for line in result.stdout.splitlines() if line.startswith("TORQUE ")]
        self.assertEqual(len(native), len(SURFACES), result.stdout)
        selected = DISSIPATIVE if ifowp or ifowt else HERMITIAN
        for cell, row in enumerate(native):
            factor = polarity * CELL_SCALE[cell]
            self.assertEqual(int(row[0]), cell + 1)
            self.assertEqual(row[1], factor * selected[0], "original ion slot")
            self.assertEqual(row[2], factor * selected[1], "original electron slot")
            accepted = selected if third_f0 == 0 else selected[:2]
            self.assertEqual(row[3], factor * sum(accepted), "all thermal species total")

        forms = rows(self.work / "NTV_TORQUE_FORMS.OUT")
        self.assertIsNotNone(forms, "original torque forms not written")
        self.assertEqual(len(forms), len(SURFACES))
        for cell, row in enumerate(forms):
            self.assertEqual(len(row), 8)
            self.assertEqual(row[:2], [cell + 1, SURFACES[cell]])
            factor = polarity * CELL_SCALE[cell]
            expected = [factor * values[species] for species in (0, 1)
                        for values in (HERMITIAN, DISSIPATIVE, WEAK)]
            self.assertEqual(row[2:], expected, "original first-two forms")

        species_rows = rows(self.work / "NTV_TORQUE_SPECIES.OUT")
        self.assertIsNotNone(species_rows, "per-species torque report not written")
        self.assertEqual(len(species_rows), 3 * len(SURFACES))
        observed = {}
        for row in species_rows:
            self.assertEqual(len(row), 7, "KP IS CSM TH TD TW SELECTED columns")
            kp, cell = int(row[0]), int(row[1])
            self.assertIn(kp, (1, 2, 3))
            self.assertIn(cell, (1, 2))
            self.assertNotIn((kp, cell), observed, "duplicate species/cell")
            observed[kp, cell] = row
            self.assertEqual(row[2], SURFACES[cell - 1])
            factor = polarity * CELL_SCALE[cell - 1]
            self.assertEqual(row[3:6], [factor * values[kp - 1]
                                      for values in (HERMITIAN, DISSIPATIVE, WEAK)])
            expected = 0.0 if kp == 3 and third_f0 != 0 else factor * selected[kp - 1]
            self.assertEqual(row[6], expected, "selected thermal contribution")
        for cell, row in enumerate(native, 1):
            self.assertEqual(sum(observed[kp, cell][6] for kp in (1, 2, 3)), row[3])

    def test_zero_width_hermitian_torque_includes_third_maxwellian(self):
        for polarity in (-1, 1):
            with self.subTest(polarity=polarity):
                self.run_case(0, 0, polarity=polarity)

    def test_each_finite_width_flag_selects_dissipative_torque(self):
        for flags in ((1, 0), (0, 1), (1, 1)):
            with self.subTest(flags=flags):
                self.run_case(*flags)

    def test_nonmaxwellian_third_slot_is_excluded_only_from_thermal_total(self):
        for flags in ((0, 0), (1, 0), (0, 1)):
            with self.subTest(flags=flags):
                self.run_case(*flags, third_f0=1)


if __name__ == "__main__":
    unittest.main()
