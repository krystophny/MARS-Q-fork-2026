"""Vacuum solve on a grid that reaches the magnetic axis, against BS_ESC.

BS_ESC (KKF=-3) fits equivalent surface currents on IFEED to a prescribed normal
field b1_m on the vacuum surface IBNM with its own Biot-Savart kernel. Driving the
coil-inductance solve (INCFEED=4, KKF=-2) with those currents must give back b1_m
on IBNM. Runs CHEASE on EXAMPLE/Data/EXPEQ_EXAMPLE at reduced resolution and MARS
on the EXAMPLE/RmpLin vacuum deck (MARS_EXE); skipped without both executables.
"""
import os
import re
import shutil
import subprocess
import unittest

from support import EXE as CHEASE, NativeTest, ROOT
from test_axis_inductance_rows import heredoc, mars, setting

NPSI = 48
MMAX = 8
KV_BNM = 10
BNM = {2: 1.0e-4, 3: -0.5e-4j}


def b1_on_surface(path, row):
    lines = path.read_text().split('\n')
    nm, ns = (int(float(x)) for x in lines[0].split()[:2])
    modes = [int(round(float(lines[1 + i].split()[0]))) for i in range(nm)]
    body = [ln.split() for ln in lines[1 + nm:] if ln.strip()]
    return {m: complex(float(body[i * ns + row][0]), float(body[i * ns + row][1]))
            for i, m in enumerate(modes)}


class AxisVacuumRoundtrip(NativeTest):
    def setUp(self):
        super().setUp()
        if mars() is None:
            raise unittest.SkipTest('set MARS_EXE for the full-program axis test')
        probe = subprocess.run([str(CHEASE)], input='', cwd=self.work, capture_output=True,
                               text=True, timeout=60) if CHEASE.exists() else None
        if probe is None or probe.returncode == 127:
            raise unittest.SkipTest('set CHEASE_EXECUTABLE to a runnable native build')

    def mars_stage(self, name, mmax, **keys):
        work = self.work / name
        work.mkdir()
        for f in ('OUTRMAR', 'OUTVMAR'):
            shutil.copy2(self.work / 'chease' / f, work / f)
        for source, target in (('PROFDEN_EXAMPLE.IN', 'PROFDEN.IN'),
                               ('PROFROT_EXAMPLE.IN', 'PROFROT.IN'),
                               ('PANGLE_EXAMPLE_PEST', 'PANGLE_PEST.IN'),
                               ('PANGLE_EXAMPLE_GEOM', 'PANGLE_EQAC.IN')):
            (work / target).write_text((ROOT / 'EXAMPLE/Data' / source).read_text())
        run = heredoc(ROOT / 'EXAMPLE/RmpLin/Rmars_EXAMPLE_RmpLin_vac')
        run = re.sub(r'&KINETIC.*?&END', '&KINETIC\n INCKIN = 0,\n&END', run, flags=re.S)
        run = re.sub(r'FEEDI\s*=[^\n]*', 'FEEDI = (1.0,0.0),', run)
        base = dict(M1=-mmax, M2=mmax, NV=40, IWALL=30, NCOIL=1, IFEED=KV_BNM + 2,
                    ISENS=KV_BNM + 2)
        for key, value in {**base, **keys}.items():
            run = setting(run, key, value)
        (work / 'RUN.IN').write_text(run + '\n')
        return work

    def run_mars(self, work):
        env = {**os.environ, 'OMP_NUM_THREADS': '1'}
        result = subprocess.run([str(mars())], cwd=work, capture_output=True, text=True,
                                timeout=900, env=env)
        output = result.stdout + result.stderr
        self.assertFalse('ZGETRF ERROR' in output, output[-1500:])
        self.assertEqual(result.returncode, 0, output[-1500:])
        return output

    def roundtrip_misfit(self, mmax):
        chease = self.work / 'chease'
        chease.mkdir()
        deck = heredoc(ROOT / 'EXAMPLE/Equi/Rchease_EXAMPLE_base')
        for key, value in dict(NS=24, NT=48, NPSI=NPSI, NCHI=128, NV=40, NPLOT=0).items():
            deck = setting(deck, key, value)
        (chease / 'EXPEQ').write_text((ROOT / 'EXAMPLE/Data/EXPEQ_EXAMPLE').read_text())
        result = subprocess.run([str(CHEASE)], input=deck + '\n', cwd=chease,
                                capture_output=True, text=True, timeout=300,
                                env={**os.environ, 'OMP_NUM_THREADS': '1'})
        self.assertEqual(result.returncode, 0, result.stdout[-2000:] + result.stderr)

        ibnm = NPSI + KV_BNM
        bnm = {m: BNM.get(m, 0j) for m in range(-mmax, mmax + 1)}
        esc = self.mars_stage('esc', mmax, KKF=-3)
        (esc / 'INPUT_BNM.IN').write_text(
            f'{len(bnm)} {ibnm}\n' + ''.join(f'{v.real: .16e} {v.imag: .16e}\n'
                                             for v in bnm.values()))
        self.run_mars(esc)
        vac = self.mars_stage('vac', mmax, KKF=-2)
        shutil.copy2(esc / 'CURHARMO.OUT', vac / 'CURHARMO.IN')
        self.run_mars(vac)
        got = b1_on_surface(vac / 'BPLASMA.OUT', ibnm - 1)
        err = sum(abs(got[m] - bnm[m]) ** 2 for m in bnm) ** 0.5
        return err / sum(abs(v) ** 2 for v in bnm.values()) ** 0.5

    def test_vacuum_solve_reproduces_the_fitted_normal_field(self):
        # At this resolution the upstream axis cell leaves a discretisation misfit of
        # 0.0102; half-cell-1 rows without curl b = 0 (fork #103 with #105) give 0.0144.
        self.assertLess(self.roundtrip_misfit(MMAX), 0.012)


if __name__ == '__main__':
    unittest.main()
