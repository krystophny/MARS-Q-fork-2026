"""Exercise the compiled vacuum limit and reject unsupported input as failure."""
import math
import os
import re
import resource
import shutil
import subprocess
import unittest

from support import EXE, ROOT, NativeTest


class VacuumCapacity(NativeTest):
    def run_deck(self, nv):
        if not EXE.exists():
            self.skipTest('set CHEASE_EXECUTABLE to a native Make build')
        deck = (ROOT/'EXAMPLE/Equi/Rchease_EXAMPLE_base').read_text()
        deck = deck.split("<<'EOF'\n", 1)[1].split('\nEOF', 1)[0]
        for key, value in dict(NS=24, NT=48, NPSI=64, NCHI=64, NV=nv,
                               NPLOT=0, NIDEAL=0).items():
            pattern = r'(?i)\b'+key+r'\s*=\s*[^,\n]+'
            if re.search(pattern, deck):
                deck = re.sub(pattern, lambda m: key+'='+str(value), deck)
            else:
                deck = deck.replace('&END', f' {key}={value},\n &END')
        shutil.copyfile(ROOT/'EXAMPLE/Data/EXPEQ_EXAMPLE', self.work/'EXPEQ')
        resource.setrlimit(resource.RLIMIT_STACK,
                           (resource.RLIM_INFINITY, resource.RLIM_INFINITY))
        return subprocess.run([str(EXE)], input=deck+'\n', cwd=self.work,
            capture_output=True, text=True, timeout=180,
            env={**os.environ, 'OMP_NUM_THREADS': '1'})

    def test_supported_limit_produces_complete_finite_vacuum_map(self):
        result = self.run_deck(640)
        self.assertEqual(result.returncode, 0, (result.stdout+result.stderr)[-2500:])
        self.assertIn('FINAL OUTPUT', result.stdout)
        self.assertNotRegex(result.stdout, r'NOT CONVERGED|NO CONVERGENCE|\bNaN\b')
        values = (self.work/'OUTVMAR').read_text().split()
        # OUTVMAR holds NV+1 radial nodes, 2*NCHI poloidal points, and 3 walls.
        self.assertEqual(list(map(int, values[:3])), [641, 128, 3])
        numbers = [float(word.replace('D', 'E')) for word in values[3:]]
        self.assertTrue(all(map(math.isfinite, numbers)))
        # Every half-interval contains one radial coordinate and seven fields.
        self.assertGreaterEqual(len(numbers), (2*640+1)*(1+7*128))

    def test_unsupported_limit_returns_nonzero_before_mapping(self):
        result = self.run_deck(641)
        self.assertGreater(result.returncode, 0, result.stderr)
        self.assertRegex(result.stdout, r'NV\s+\.GT\.\s+NPV')
        self.assertNotIn('FINAL OUTPUT', result.stdout)
        self.assertFalse((self.work/'OUTVMAR').exists())


if __name__ == '__main__':
    unittest.main()
