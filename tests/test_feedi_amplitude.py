"""KFEEDIAMP fixes what FEEDI means; every torque is then the physical one.

For a real coil current I0 cos(n phi), FEEDI = I0 / c with c = KFEEDIAMP:
c = 1 is the cos-sin (peak) coefficient, c = 2 the +n coefficient of the complex Fourier
series. MARS fields are per FEEDI, so the physical field is c Re(b exp(i n phi)) and every
quadratic output carries c^2 on top of its peak-amplitude formula (1/2 Re(J* B) for j x B,
-RNTOR Im(dW_K) with dW_K = 1/2 int xi* F for NTV).

PhysicsOracle always runs. KernelScaling needs MARS_EXE and MARS_RMP_CASE, a directory with
a complete MARS-F RMP run that computes the j x B, NTV and Reynolds torques (for example
EXAMPLE/RmpLin pls with its OUT*MAR, PROF*.IN and PANGLE*.IN files and RUN.IN). It runs the
case with KFEEDIAMP = 1 and 2: the fields must be identical and every total torque must
scale by exactly 4.
"""

from __future__ import annotations

import cmath
import math
import os
import re
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

TORQUES = ("JXB", "NTV", "REYNOLDS")


def _phase_average(fn, nphi: int = 4096) -> float:
    return sum(fn(2.0 * math.pi * k / nphi) for k in range(nphi)) / nphi


class PhysicsOracle(unittest.TestCase):
    def test_torque_of_real_coil_current_for_both_conventions(self) -> None:
        n, i0 = -3, 1.7  # n plays RNTOR: fields vary as exp(i RNTOR phi)
        hhat, fhat = complex(0.7, -0.2), complex(-0.3, 0.9)  # response per unit FEEDI
        physical = None
        for c in (1, 2):
            feedi = i0 / c
            h, f = feedi * hhat, feedi * fhat
            direct = _phase_average(lambda p: -(c * f * cmath.exp(1j * n * p)).real
                                    * (c * 1j * n * h * cmath.exp(1j * n * p)).real)
            dwk = 0.5 * h.conjugate() * f
            self.assertAlmostEqual(direct, -n * c**2 * dwk.imag, places=12)
            self.assertNotAlmostEqual(direct, -2 * n * dwk.imag, places=3)  # old coefficient
            if physical is None:
                physical = direct
            self.assertAlmostEqual(direct, physical, places=12)  # same coil current

    def test_jxb_form_is_peak_form_times_c_squared(self) -> None:
        j, b = complex(0.4, 1.1), complex(-0.8, 0.3)
        for c in (1, 2):
            direct = _phase_average(lambda p: (c * j * cmath.exp(1j * p)).real
                                    * (c * b * cmath.exp(1j * p)).real)
            self.assertAlmostEqual(direct, c**2 * 0.5 * (j.conjugate() * b).real, places=12)


@unittest.skipUnless(os.environ.get("MARS_EXE") and os.environ.get("MARS_RMP_CASE"),
                     "set MARS_EXE and MARS_RMP_CASE for the KFEEDIAMP runs")
class KernelScaling(unittest.TestCase):
    def test_torques_scale_with_kfeediamp_squared(self) -> None:
        exe = Path(os.environ["MARS_EXE"]).resolve()
        src = Path(os.environ["MARS_RMP_CASE"])
        runs = {c: _run(exe, src, c) for c in (1, 2)}
        self.assertEqual(runs[1][1], runs[2][1], "BPLASMA.OUT must not depend on KFEEDIAMP")
        for name in TORQUES:
            t1, t2 = (_totals(runs[c][0], name) for c in (1, 2))
            self.assertTrue(t1, name)
            self.assertEqual(len(t1), len(t2), name)
            for a, b in zip(t1, t2):
                self.assertNotEqual(a, 0.0, name)
                self.assertAlmostEqual(b / a, 4.0, delta=1e-4, msg=name)  # 6 printed digits


def _totals(log: str, name: str) -> list[float]:
    return [float(v) for v in re.findall(rf"TOTAL {name} TORQUE:\s*(\S+)", log)]


def _run(exe: Path, src: Path, c: int) -> tuple[str, bytes]:
    with tempfile.TemporaryDirectory(prefix=f"mars-feediamp-{c}-") as tmp:
        work = Path(tmp)
        for item in src.iterdir():
            if item.is_file():
                shutil.copy(item, work / item.name)
        run_in = (work / "RUN.IN").read_text()
        run_in = re.sub(r"\s*KFEEDIAMP\s*=\s*\d+\s*,", "", run_in)
        run_in = re.sub(r"(&FEEDBACK\s*\n)", rf"\1 KFEEDIAMP = {c},\n", run_in, count=1)
        (work / "RUN.IN").write_text(run_in)
        done = subprocess.run([str(exe)], cwd=work, stdin=subprocess.DEVNULL, capture_output=True,
                              text=True, timeout=3600, check=True)
        return done.stdout, (work / "BPLASMA.OUT").read_bytes()


if __name__ == "__main__":
    unittest.main()
