"""Reset of the accumulated kinetic surface coefficients (MPIENV::INIT_SURFACE_QUANTITY).

KJPCOEFF adds each surface contribution to the V* coefficient arrays, so every
new sweep (KJPKEY=2) and every MPI task of a worker must start from zero in all
of them: the fluid-like parallel/perpendicular blocks, the pressure-drift
(VDP*) blocks, the perturbed-potential (*DPHI*) blocks and the kinetic-electron
(*PARE/*PERE) blocks.  The production routine is compiled with the owning
module and must leave no array nonzero.
"""
import re

from support import ROOT, NativeTest

NATIVE = ROOT / "MarsQ_2FK"


def init_surface_quantity():
    match = re.search(r"^ +SUBROUTINE +INIT_SURFACE_QUANTITY\b.*?^ +END SUBROUTINE +INIT_SURFACE_QUANTITY\b",
                      (NATIVE / "mpienv.f").read_text(), re.M | re.S | re.I)
    return match.group(0) + "\n"


def names():
    """Every complex surface-coefficient array declared in KINETICM."""
    source = (NATIVE / "kineticm.f").read_text()
    blocks = re.findall(r"COMPLEX\*16,DIMENSION\(:,:,:\),ALLOCATABLE::(.*?)\n(?! +&)",
                        source, re.S)
    found = []
    for block in blocks[:2]:
        found += [n for n in re.split(r"[\s,&]+", block) if n]
    return found


class SurfaceReset(NativeTest):
    def test_new_sweep_starts_from_zero_in_every_coefficient(self):
        arrays = names()
        self.assertIn("VDPPARA", arrays)
        self.assertIn("VX1DPHI", arrays)
        self.assertIn("VDPPEREM", arrays)
        lines = ["      SUBROUTINE CHECK_RESET", "      USE KINETICM", "      IMPLICIT NONE",
                 "      INTEGER BAD", "      BAD = 0"]
        for name in arrays:
            lines += ["      ALLOCATE(%s(3,3,2))" % name, "      %s = (1.0,2.0)" % name]
        lines.append("      CALL INIT_SURFACE_QUANTITY()")
        for name in arrays:
            lines.append("      IF (ANY(%s.NE.(0.0,0.0))) THEN" % name)
            lines.append("         WRITE(*,*) 'NOT RESET: %s'" % name)
            lines.append("         BAD = BAD + 1")
            lines.append("      ENDIF")
        lines += ["      IF (BAD.GT.0) STOP 1", "      END SUBROUTINE CHECK_RESET",
                  "      PROGRAM SURFACE_RESET", "      CALL CHECK_RESET", "      END PROGRAM"]
        exe = self.compile({
            "globalm.f": (NATIVE / "globalm.f").read_text(),
            "kineticm.f": (NATIVE / "kineticm.f").read_text(),
            "reset.f": init_surface_quantity(),
            "main.f": "\n".join(lines) + "\n",
        })
        result = self.execute(exe)
        self.assertEqual(result.stdout.strip(), "")


if __name__ == "__main__":
    import unittest
    unittest.main()
