"""PROFWE = 0 must give a finite NTV torque (no 0/0 when normalising the profile).

The source contract runs everywhere.  The executable regression needs MARS_EXE
and MARS_NTV_CASE, a directory with a complete Shaing-NTV run (RUN.IN with
KNTV=10, KEYTORQ=1, NPROFWE=4, INCKIN=0; PROF*.IN, PROFWE.IN, OUTRMAR, OUTVMAR);
the test zeroes PROFWE.IN and runs MARS in a copy.
"""

from __future__ import annotations

import math
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "MarsQ_2FK"


class SourceContract(unittest.TestCase):
    def test_profile_normalisation_guards_zero_first_value(self) -> None:
        text = (SRC / "marsq.f").read_text(errors="replace")
        self.assertRegex(text, r"IF \(ZTEMP\.EQ\.0\.0\) THEN\s+IF \(MAXVAL\(ABS\(RPRV\)\)\.GT\.0\.0\) STOP")
        self.assertNotRegex(text, r"RPRV\s+=\s+ROTWE0\*RPRV/ZTEMP\s*\n\s+IF \(IRAD")


@unittest.skipUnless(os.environ.get("MARS_EXE") and os.environ.get("MARS_NTV_CASE"),
                     "set MARS_EXE and MARS_NTV_CASE for the PROFWE=0 run")
class ZeroRotation(unittest.TestCase):
    def test_torque_is_finite(self) -> None:
        exe = Path(os.environ["MARS_EXE"]).resolve()
        src = Path(os.environ["MARS_NTV_CASE"])
        with tempfile.TemporaryDirectory(prefix="mars-profwe0-") as tmp:
            work = Path(tmp)
            for item in src.iterdir():
                if item.is_file():
                    shutil.copy(item, work / item.name)
            lines = (work / "PROFWE.IN").read_text().splitlines()
            rows = [f"{r.split()[0]} 0.0" for r in lines[1:] if r.strip()]
            (work / "PROFWE.IN").write_text("\n".join([lines[0], *rows]) + "\n")
            subprocess.run([str(exe)], cwd=work, stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                           stderr=subprocess.STDOUT, timeout=1800, check=True)
            values = [float(v) for line in (work / "TORQUENTV.OUT").read_text().splitlines()
                      for v in line.split()]
            self.assertTrue(values and all(math.isfinite(v) for v in values))
            self.assertTrue(any(v != 0.0 for v in values))


if __name__ == "__main__":
    unittest.main()
