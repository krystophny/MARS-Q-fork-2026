"""Opt-in wall reload preserves the native plasma and rejects changed contours."""
import math
import subprocess
import os
import re
import shutil
import struct
import unittest
from support import NativeTest, NATIVE, EXE, ROOT, routine


class RestartWall(NativeTest):
    def setUp(self):
        super().setUp()
        files = {p.name: p.read_text() for p in NATIVE.glob('*.inc')}
        files['native.f'] = '\n'.join(routine(n) for n in
            ('RSTWALL', 'BNDSPLW', 'SPLCY', 'TRICYC', 'SORT',
             'SCOPY', 'ISMAX', 'ISMIN'))
        files['oracle.f'] = '''
      PROGRAM ORACLE
      INCLUDE 'DECLAR.inc'
      INCLUDE 'COMDIM.inc'
      INCLUDE 'COMBND.inc'
      INCLUDE 'COMCON.inc'
      INCLUDE 'COMESH.inc'
      INCLUDE 'COMSOL.inc'
      DIMENSION ZSAVE(17,5),ZWALL(17,2)
      NINWALL=1
      NOPT=1
      NSURF=6
      NSYM=0
      NBPS=17
      NWBPS=2
      RZ0=0.07
      R0W=0.9
      RZ0W=0.
      REXT=5.
      BPS(1)=3.
      BPS(12)=0.
      CPI=ACOS(-1.)
      RC2PI=2.*CPI
      RC1M14=1.D-14
      DO L=1,NBPS
         T=2.*CPI*(L-1)/(NBPS-1)
         RRBPS(L,1)=3.+0.8*COS(T)
         RZBPS(L,1)=0.8*SIN(T)
         TETBPS(L,1)=-L
         D2RBPS(L,1)=17.*L
         D2ZBPS(L,1)=-23.*L
      ENDDO
      ZSAVE(:,1)=RRBPS(1:NBPS,1)
      ZSAVE(:,2)=RZBPS(1:NBPS,1)
      ZSAVE(:,3)=TETBPS(1:NBPS,1)
      ZSAVE(:,4)=D2RBPS(1:NBPS,1)
      ZSAVE(:,5)=D2ZBPS(1:NBPS,1)
      CALL RSTWALL
      IF (NBPS.NE.17.OR.NWBPS.NE.3) STOP 2
      IF (ANY(ZSAVE(:,1).NE.RRBPS(1:NBPS,1))) STOP 3
      IF (ANY(ZSAVE(:,2).NE.RZBPS(1:NBPS,1))) STOP 4
      IF (ANY(ZSAVE(:,3).NE.TETBPS(1:NBPS,1))) STOP 5
      IF (ANY(ZSAVE(:,4).NE.D2RBPS(1:NBPS,1))) STOP 6
      IF (ANY(ZSAVE(:,5).NE.D2ZBPS(1:NBPS,1))) STOP 7
      ERR=ABS(RW(2)-(4.2-0.9)/(3.8-0.9))
      ERR=MAX(ERR,ABS(RW(3)-5.))
      DO L=1,NBPS
         T=TETBPS(L,2)
         ERR=MAX(ERR,ABS(RRBPS(L,2)-3.-1.2*COS(T)),
     &               ABS(RZBPS(L,2)-1.2*SIN(T)))
      ENDDO
      PRINT *, 'wall analytic error=',ERR
      IF(ERR.GT.1.D-12) STOP 8
      IF (ANY(CNDRZ(1:NBPS,2:NWBPS).NE.1.)) STOP 9
      IF (MAXVAL(ABS(D2CBPS(1:NBPS,2:NWBPS))).GT.1.D-12) STOP 10
      PRINT *, 'native plasma arrays unchanged'
      END
'''
        self.exe = self.compile(files)

    def write_geometry(self, shift=0., count=17, radius=1.2):
        rows = ['0.2666666666666667', '0.07', '0.', f'{count} 2 2']
        for rad, delta in ((.8, shift), (radius, 0.)):
            for i in range(17):
                angle = 2*math.pi*i/16
                rows.append(f'{3+rad*math.cos(angle)+delta:.17e} {rad*math.sin(angle):.17e}')
        (self.work/'EXPEQ').write_text('\n'.join(rows)+'\n')

    def test_wall_reload_preserves_plasma_and_matches_analytic_circle(self):
        self.write_geometry()
        self.assertIn('native plasma arrays unchanged', self.execute(self.exe).stdout)

    def test_changed_plasma_contour_rejected(self):
        self.write_geometry(shift=.005)
        result = subprocess.run([str(self.exe)], cwd=self.work, capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('PLASMA CONTOUR DIFFERS', result.stdout)

    def test_inconsistent_vertex_count_rejected(self):
        self.write_geometry(count=16)
        result = subprocess.run([str(self.exe)], cwd=self.work, capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('INVALID EXPEQ WALL GEOMETRY', result.stdout)

    def test_non_enclosing_wall_rejected(self):
        self.write_geometry(radius=.7)
        result = subprocess.run([str(self.exe)], cwd=self.work, capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('INVALID EXPEQ WALL GEOMETRY', result.stdout)


class RestartWallFullProgram(NativeTest):
    @unittest.skipUnless(EXE.exists(), 'set CHEASE_EXECUTABLE to a native Make build')
    def test_frozen_native_state_and_mapping_match_existing_wall_route(self):
        lines = (ROOT/'EXAMPLE/Data/EXPEQ_EXAMPLE').read_text().splitlines()
        nb = int(lines[3].split()[0])
        full = '\n'.join(lines)+'\n'
        lines[3] = f'{nb} 1 2'
        del lines[4+nb:4+2*nb]
        single = '\n'.join(lines)+'\n'
        deck = (ROOT/'EXAMPLE/Equi/Rchease_EXAMPLE_base').read_text()
        deck = deck.split("<<'EOF'\n", 1)[1].split('\nEOF', 1)[0]

        def setting(text, key, value):
            if re.search(r'(?i)\b'+key+r'\s*=', text):
                return re.sub(r'(?i)\b'+key+r'\s*=\s*[^,\n]+',
                              lambda m: key+'='+str(value), text)
            return text.replace('&END', f' {key}={value},\n &END')

        for key, value in dict(NS=24, NT=48, NPSI=64, NCHI=64, NV=40,
                               NPLOT=0, NIDEAL=6).items():
            deck = setting(deck, key, value)

        def execute(label, text, geometry, checkpoint=None):
            work = self.work/label
            work.mkdir()
            (work/'EXPEQ').write_text(geometry)
            if checkpoint:
                shutil.copyfile(checkpoint, work/'NIN')
            result = subprocess.run([str(EXE)], input=text+'\n', cwd=work,
                capture_output=True, text=True, timeout=120,
                env={**os.environ, 'OMP_NUM_THREADS': '1'})
            self.assertEqual(result.returncode, 0, result.stdout[-2500:]+result.stderr)
            self.assertNotRegex(result.stdout, r'NOT CONVERGED|NO CONVERGENCE|\bNaN\b')
            if checkpoint:
                self.assertEqual((work/'NIN').read_bytes(), checkpoint.read_bytes())
                self.assertNotRegex(result.stdout, r'RESIDU\s*=')
            return work

        one = execute('single', deck, single)
        two = execute('with_wall', deck, full)

        def records(path):
            raw, offset, out = path.read_bytes(), 0, []
            while offset < len(raw):
                count = struct.unpack_from('<i', raw, offset)[0]
                self.assertEqual(struct.unpack_from('<i', raw, offset+count+4)[0], count)
                out.append(raw[offset+4:offset+4+count])
                offset += count+8
            return out

        before, after = records(one/'NOUT'), records(two/'NOUT')
        self.assertEqual(len(before), 47)
        self.assertEqual(len(after), 47)
        for index in range(47):
            if not 34 <= index <= 40:
                self.assertEqual(before[index], after[index], f'physical record {index}')
        for index in range(36, 41):
            self.assertEqual(before[index][:8*nb], after[index][:8*nb],
                             f'native plasma boundary record {index}')
        mapping = setting(setting(deck, 'NOPT', 1), 'NIDEAL', 0)
        reference = execute('native_wall', mapping, full, two/'NOUT')
        candidate = execute('wall_reload', setting(mapping, 'NINWALL', 1), full, one/'NOUT')
        for filename in ('OUTRMAR', 'OUTVMAR'):
            self.assertEqual((reference/filename).read_bytes(),
                             (candidate/filename).read_bytes(), filename)
        print('47-record native physical-state parity; frozen OUTRMAR/OUTVMAR byte equality')
