"""KNTV = 20/21 torque is T = n Im(dW_K) (dW_K already carries the 1/2), not 2 n Im(dW_K).

Convention: perturbations are Re(X exp(i n phi)) (torque.f takes Re(J^* B)/2) and the MARS-K
energy is dW = 1/2 int H^* f, so the torque, T = -<f dH/dphi> averaged over the real fields,
is n Im(dW) (per unit (2 pi)^2 of the chi, phi angle integration).

PhysicsOracle always runs.  KernelTorque needs MARS_EXE and MARS_NTV_CASE, a directory with a
complete native perturbative MARS-K run (RUN.IN with INCKIN = 1, IPERTURB = 2, ODWKCOM =
.false., KNTV = 20, a perturbation drive such as INCFEED = 8, PROF*.IN, OUTRMAR, OUTVMAR).
It runs the case twice in copies: KNTV = 20 prints "TOTAL NTV TORQUE" (the radial sum of
TORQUENTV, 4 pi^2 sum TORQUENTV CSH) and KNTV = 21 prints "TOTAL KINETIC ENERGY" (dW_K of
the same response); the torque must equal -RNTOR Im(dW_K), with n = -RNTOR.
KNTV = 21 itself needs the DATA_STORAGE/PCOEF files of the fixed-field path (issues 2, 3) and
cannot run on a native response.
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


class PhysicsOracle(unittest.TestCase):
    def test_real_field_phase_average_gives_n_im_dw(self) -> None:
        n, nphi = 3, 4096
        h, f = complex(0.7, -0.2), complex(-0.3, 0.9)  # H = Re(h e^{i n phi}), f = chi h
        dw = 0.5 * h.conjugate() * f
        tsum = 0.0
        for k in range(nphi):
            phi = 2.0 * math.pi * k / nphi
            fr = (f * cmath.exp(1j * n * phi)).real
            dh = (1j * n * h * cmath.exp(1j * n * phi)).real
            tsum += -fr * dh
        direct = tsum / nphi
        self.assertAlmostEqual(abs(direct), n * abs(dw.imag), places=12)
        self.assertGreater(abs(2.0 * n * dw.imag), 1.5 * abs(direct))  # the pre-fix coefficient


@unittest.skipUnless(os.environ.get("MARS_EXE") and os.environ.get("MARS_NTV_CASE"),
                     "set MARS_EXE and MARS_NTV_CASE for the KNTV=20/21 runs")
class KernelTorque(unittest.TestCase):
    def test_total_torque_is_n_im_dwk(self) -> None:
        exe = Path(os.environ["MARS_EXE"]).resolve()
        src = Path(os.environ["MARS_NTV_CASE"])
        out = {kntv: _run(exe, src, kntv) for kntv in (20, 21)}
        rntor = float(re.search(r"RNTOR\s*=\s*([-+0-9.eEdD]+)", (src / "RUN.IN").read_text()).group(1))
        torque = float(re.search(r"TOTAL NTV TORQUE:\s*(\S+)", out[20]).group(1))
        dwk = re.search(r"TOTAL KINETIC ENERGY =\s*\(\s*(\S+?),\s*(\S+?)\)", out[21])
        expected = -rntor * float(dwk.group(2))
        self.assertNotEqual(expected, 0.0)
        self.assertAlmostEqual(torque / expected, 1.0, delta=2e-3)


def _run(exe: Path, src: Path, kntv: int) -> str:
    with tempfile.TemporaryDirectory(prefix=f"mars-dwk-torque-{kntv}-") as tmp:
        work = Path(tmp)
        for item in src.iterdir():
            if item.is_file():
                shutil.copy(item, work / item.name)
        run_in = re.sub(r"KNTV\s*=\s*\d+", f"KNTV = {kntv}", (work / "RUN.IN").read_text())
        (work / "RUN.IN").write_text(run_in)
        done = subprocess.run([str(exe)], cwd=work, stdin=subprocess.DEVNULL, capture_output=True,
                              text=True, timeout=3600, check=True)
        return done.stdout


if __name__ == "__main__":
    unittest.main()
