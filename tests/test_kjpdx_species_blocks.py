"""Ion and kinetic-electron blocks of the matrix-free kinetic product (KJPDX).

KJP assembles the kinetic pressure rows with the same couplings for both
species: the ion coefficients V*PARA/V*PERP fill rows KYPPARA/KYPPERP, the
electron coefficients V*PARE/V*PERE fill rows KYPE/KYPP, each pressure-drift
coefficient VDP* in the row of its own species and component.  KJPDX is the
same operator applied to a vector (IPERTURB=0).  With electron coefficients
equal to the ion ones the electron rows must therefore equal the ion rows.
"""
import re

from support import ROOT, NativeTest

NATIVE = ROOT / "MarsQ_2FK"


def kjpdx():
    match = re.search(r"^ +SUBROUTINE +KJPDX\b.*?^ +END *$", (NATIVE / "kinetic.f").read_text(),
                      re.M | re.S | re.I)
    return match.group(0) + "\n"


ION = ["X1", "X2", "Q1", "Q2", "Q3", "DP"]

DRIVER = """      PROGRAM KJPDX_BLOCKS
      USE DIMENSIM
      USE GLOBALM
      USE RCOMDM
      USE KINETICM
      IMPLICIT NONE
      INTEGER, PARAMETER :: MD = 16, MDY = 16, ND = 3
      COMPLEX*16 R(MD,ND+1), X(MD,ND+1), RY(MDY,ND), Y(MDY,ND)
      REAL*8 A(MDY,ND,2), B(MD,ND+1,2)
      INTEGER I, M
      NR = 3
      MSMAX = 2
      NXCOMP = 2
      NYCOMP = 8
      KXB1 = 1
      KXDPHI = 2
      KYB2 = 1
      KYB3 = 2
      KYPPARA = 3
      KYPPERP = 4
      KYPE = 5
      KYPP = 6
      IDIAMTE = 0
      IDIAMB = 1
      INCKIN = 1
      PTRAPH = 0.9
      CALPHA1 = (1.3,0.2)
      CALPHA2 = (0.7,-0.1)
      CALPHA4 = (0.4,0.5)
      ALLOCATE(CS(NR+1), CSM(NR), CSH(NR))
      CS = [0.1, 0.4, 0.7, 1.0]
      CSM = 0.5*(CS(1:NR)+CS(2:NR+1))
      CSH = CS(2:NR+1)-CS(1:NR)
{allocate}
      R = 0
      RY = 0
      CALL RANDOM_SEED(PUT=[(7*I+3, I=1,64)])
      CALL RANDOM_NUMBER(A)
      CALL RANDOM_NUMBER(B)
      Y = DCMPLX(A(:,:,1),A(:,:,2))
      X = DCMPLX(B(:,:,1),B(:,:,2))
      CALL KJPDX(MD, MDY, ND, R, RY, X, Y)
      DO I = 1, ND
         DO M = 0, MSMAX-1
            WRITE(*,'(8ES24.15)') RY(KYPPARA+M*NYCOMP,I),
     &         RY(KYPE+M*NYCOMP,I), RY(KYPPERP+M*NYCOMP,I),
     &         RY(KYPP+M*NYCOMP,I)
         END DO
      END DO
      END PROGRAM
"""


def allocate():
    lines = ["      CALL RANDOM_SEED(PUT=[(5*I+1, I=1,64)])"]
    for name in ION:
        for comp, ecomp in (("PARA", "PARE"), ("PERP", "PERE")):
            for full, mesh in ((name + comp, "NR+1"), (name + comp + "M", "NR")):
                efull = "V" + full.replace(comp, ecomp)
                full = "V" + full
                lines += ["      ALLOCATE(%s(MSMAX,MSMAX,%s), %s(MSMAX,MSMAX,%s))"
                          % (full, mesh, efull, mesh),
                          "      CALL RANDOM_NUMBER(A(1:MSMAX*MSMAX*(%s),1,1))" % mesh,
                          "      %s = RESHAPE(DCMPLX(A(1:MSMAX*MSMAX*(%s),1,1),0.3D0),"
                          " SHAPE(%s))" % (full, mesh, full),
                          "      %s = %s" % (efull, full)]
    return "\n".join(lines)


class KjpdxSpeciesBlocks(NativeTest):
    def test_electron_rows_equal_ion_rows_for_equal_coefficients(self):
        exe = self.compile({
            "dimensim.f": (NATIVE / "dimensim.f").read_text(),
            "globalm.f": (NATIVE / "globalm.f").read_text(),
            "rcomdm.f": (NATIVE / "rcomdm.f").read_text(),
            "kineticm.f": (NATIVE / "kineticm.f").read_text(),
            "integc.inc": (NATIVE / "integc.inc").read_text(),
            "tophat.inc": (NATIVE / "tophat.inc").read_text(),
            "kjpdx.f": kjpdx(),
            "main.f": DRIVER.replace("{allocate}", allocate()),
        })
        rows = [list(map(float, line.split()))
                for line in self.execute(exe).stdout.splitlines() if line.strip()]
        self.assertEqual(len(rows), 6)
        for row in rows:
            para, pare = complex(*row[0:2]), complex(*row[2:4])
            perp, pere = complex(*row[4:6]), complex(*row[6:8])
            self.assertGreater(abs(para), 1e-3)
            self.assertGreater(abs(perp), 1e-3)
            self.assertAlmostEqual(abs(pare - para), 0.0, delta=1e-12 * abs(para))
            self.assertAlmostEqual(abs(pere - perp), 0.0, delta=1e-12 * abs(perp))


if __name__ == "__main__":
    import unittest
    unittest.main()
