"""Compile and exercise the native energy kernels, including Maxwellian tails."""
import math
import re
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "MarsQ_2FK"
FLAGS = ["-O2", "-fdefault-real-8", "-fdefault-double-8", "-fopenmp",
         "-ffixed-line-length-none", "-ffree-line-length-none",
         "-fallow-argument-mismatch", "-fcheck=all", "-ffunction-sections"]


def build_fixture(directory, anisotropic=None):
    text = Path(anisotropic or SOURCE / "anisotropic.f").read_text()
    routines = []
    for name in ("KDISTRIBF_TYPE0", "KENERGYMESH", "KIA_TRAP", "KIA_PASS",
                 "KIA_TRAP0", "KIA_PASS0"):
        # Compile the complete native routine, without replacing its arithmetic.
        match = re.search(r"^      SUBROUTINE +" + name + r"\b", text, re.M)
        start = match.start()
        end = text.index("END SUBROUTINE " + name, start)
        routines.append(text[start:end] + "END SUBROUTINE " + name + "\n")
    kernels = Path(directory) / "kernels.f"
    kernels.write_text("\n".join(routines))
    command = ["gfortran", *FLAGS, "-I", str(SOURCE),
               *[str(SOURCE / name) for name in
                 ("dimensim.f", "globalm.f", "rcomdm.f", "kineticm.f",
                  "anisotropicm.f")], str(kernels), str(SOURCE / "splines.f"),
               str(ROOT / "tests/dwk_energy_tail.f90"),
               "-Wl,--gc-sections", "-o", "fixture"]
    subprocess.run(command, cwd=directory, check=True, capture_output=True, text=True)
    return Path(directory) / "fixture"


def run_fixture(executable, nepk):
    return subprocess.check_output([str(executable), str(nepk)], text=True)


def values(output):
    result = {}
    for line in output.splitlines():
        fields = line.split()
        if fields[0] in ("TRAP", "PASS", "TRAP0", "PASS0"):
            result[tuple(fields[:7])] = complex(float(fields[7]), float(fields[8]))
    return result


def maxwellian_oracle(nepk, route, pitch, kp, harmonic, coll):
    """Direct Gaussian sum of analytic Maxwellian derivatives and resolvent."""
    pi = 3.14159265358979
    a = -3.0*1.1/50.0*0.9 * (1 if kp == 1 else -1)
    omega0 = complex(-0.09, -1e-11)
    mass = 2.0 if kp == 1 else 0.00054463
    bounce = 0.4 + 0.1*pitch
    drift = -0.3 + 0.2*pitch
    q = 2.1
    collision = 0.05
    total = 0j
    absolute_sum = 0.0
    # Accumulate interval edges exactly as KENERGYMESH does.
    h = 1.0/(nepk-1)
    left = 0.0
    for interval in range(nepk-1):
        right = left+h
        for w in (1/math.sqrt(3), -1/math.sqrt(3)):
            u = 0.5*(left*(1+w)+right*(1-w))
            energy = math.tan(u*0.5*math.pi)
            f = 1.3*math.exp(-energy)
            fpsi = f*(-0.7+(energy-2.5)*(-0.4))
            numerator = (a*fpsi-omega0*f)*energy**2.5
            jacobian = 0.5*pi/math.cos(0.5*pi*u)**2
            scale = 1.0 if route == "TRAP" else 0.9**(-1.5)
            harmonic_scale = 1.0+harmonic**2/4.0 if coll == 1 else 1.0
            damping = complex(0, -collision*energy**(-1.5)*scale*harmonic_scale)
            if route == "TRAP":
                frequency = a*drift*energy if harmonic == 0 else (
                    bounce*harmonic*math.sqrt(2*0.9*2.0/mass)*math.sqrt(energy)
                    +a*drift*energy)
                value = -2*numerator/(omega0+frequency+damping)
                magnitude = abs(value)
            else:
                frequency = bounce*(-3*q+harmonic)*math.sqrt(2*0.9*2.0/mass)
                co = -numerator/(omega0+frequency*math.sqrt(energy)+damping)
                counter = -numerator/(omega0-frequency*math.sqrt(energy)+damping)
                value = co+counter
                magnitude = abs(co)+abs(counter)
            total += 0.5*(right-left)*jacobian*value
            absolute_sum += 0.5*(right-left)*jacobian*magnitude
        left = right
    return total, absolute_sum


@unittest.skipUnless(shutil.which("gfortran"), "gfortran required")
class EnergyTailTests(unittest.TestCase):
    def test_native_kernels_against_maxwellian_oracle(self):
        with tempfile.TemporaryDirectory() as directory:
            executable = build_fixture(directory)
            for nepk in (201, 401, 801):
                result = values(run_fixture(executable, nepk))
                self.assertEqual(len(result), 1728)
                for key, value in result.items():
                    self.assertTrue(math.isfinite(value.real) and math.isfinite(value.imag), key)
                    route, coll, grid, pitch, kp, opt, harmonic = key
                    if coll in ("1", "2") and opt == "0" and route in ("TRAP", "PASS"):
                        expected, scale = maxwellian_oracle(
                            nepk, route, int(pitch), int(kp), int(harmonic)-2, int(coll))
                        # Co/counter terms can cancel strongly for electrons.
                        self.assertLess(abs(value-expected),
                                        64*sys.float_info.epsilon*scale, key)


if __name__ == "__main__":
    unittest.main()
