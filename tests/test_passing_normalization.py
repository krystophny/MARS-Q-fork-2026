"""Native passing response must be invariant under reference-density units.

Compile complete, unchanged production KIA_PASS/KIA_PASS0 routines with the
actual owning module declarations and native spline implementation. The oracle
changes units while keeping the physical DKE fixed; it does not duplicate the
passing integral or inspect patch text.
"""
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile
import unittest

ROOT = Path(os.environ.get('MARSQ_ROOT') or Path(__file__).resolve().parents[1])
NATIVE = ROOT / 'MarsQ_2FK'


def production_routine(name):
    source = (NATIVE / 'anisotropic.f').read_text()
    match = re.search(r'^ +SUBROUTINE +' + name + r'\b.*?^ +END SUBROUTINE +' + name + r'\b',
                      source, re.MULTILINE | re.DOTALL | re.IGNORECASE)
    if match is None:
        raise ValueError('Missing production routine ' + name)
    return match.group() + '\n'


class PassingNormalizationTests(unittest.TestCase):
    def test_reference_density_does_not_change_passing_response(self):
        with tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR', ROOT)) as tmp:
            work = Path(tmp)
            routines = work / 'passing.f'
            routines.write_text(production_routine('KIA_PASS') + production_routine('KIA_PASS0') +
                                "      SUBROUTINE SET_TEST_OMEGA()\n"
                                "      IMPLICIT NONE\n"
                                "      INCLUDE 'compam.inc'\n"
                                "      OMEGA=(0.D0,0.D0)\n"
                                "      END\n")
            modules = ['globalm.f', 'dimensim.f', 'rcomdm.f', 'kineticm.f', 'anisotropicm.f']
            command = [os.environ.get('FC', 'gfortran'), '-O0', '-fcheck=all', '-std=legacy',
                       '-fdefault-real-8', '-fdefault-double-8', '-ffixed-line-length-none',
                       '-ffree-line-length-none', '-fallow-argument-mismatch', '-I' + str(NATIVE),
                       *map(str, (NATIVE / name for name in modules)), str(routines),
                       str(NATIVE / 'splines.f'),
                       str(Path(__file__).with_name('fixtures') / 'passing_normalization.f90'),
                       '-o', str(work / 'passing')]
            compiled = subprocess.run(command, cwd=work, text=True, capture_output=True, timeout=60)
            self.assertEqual(compiled.returncode, 0, compiled.stdout + compiled.stderr)
            run = subprocess.run([str(work / 'passing')], cwd=work, text=True,
                                 capture_output=True, timeout=30)
            self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
            values = {}
            for line in run.stdout.splitlines():
                mode, index, eps, real, imag = line.split()
                value = complex(float(real), float(imag))
                self.assertTrue(math.isfinite(abs(value)), line)
                values.setdefault(int(mode), []).append(value)
            self.assertEqual(set(values), {0, 1, 2, 3})
            errors = {}
            for mode, result in values.items():
                self.assertEqual(len(result), 3)
                self.assertGreater(abs(result[0]), 1.e-8, 'vacuous response')
                errors[mode] = max(abs(x-result[0])/abs(result[0]) for x in result[1:])
            print('passing normalization relative errors:', errors)
            for mode, error in errors.items():
                with self.subTest(inutype=mode if mode < 3 else 'KIA_PASS0'):
                    self.assertLess(error, 1.e-11, run.stdout)


if __name__ == '__main__':
    unittest.main()
