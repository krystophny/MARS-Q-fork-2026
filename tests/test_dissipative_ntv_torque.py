#!/usr/bin/env python3
"""Native KNTV=21 torque from the DWK component contraction.

The NTV torque is n Im(dW_k).  For a purely reactive resonance kernel it must
vanish, also when the left (displacement/strain) drive of the weak
contraction differs from the right (perturbed-field) drive, as it does for a
discrete field that does not satisfy ideal induction exactly.  When both drives
agree, the torque is the known dissipative value n Im(K) c |H|^2.

Requires gfortran.  KJPFILL, the DWK component cache writer/reader and
CALCDWKCOMP with its helpers are extracted unchanged from the selected source.
"""
import argparse
import cmath
import math
from pathlib import Path
import re
import subprocess
import tempfile
import unittest

ROUTINES = ['KJPFILL', 'SETDWKCOMPMAT', 'KDWKDISSIPATIVE', 'ALLOCATEDWKCOMPMAT',
            'DEALLOCATEDWKCOMPMAT', 'WRITE_SURFACE_QUANTITIES',
            'READ_SURFACE_QUANTITIES', 'FILLMATDWKCOMP', 'CALCPRECOMP',
            'CALCDWKPROF', 'CALCDWKCOMP', 'WRITEDWKBREAKDOWN',
            'WRITEDWKDRIVELEDGER', 'CALCDWKDIRECTCHECK']

# Manufactured single-harmonic, single-orbit data.  G: parallel pressure
# projection; H: drive coefficients of (X1, X2, B1, B2, B3).
G = 0.8 + 0.3j
H = [0.5 - 0.2j, 0.3 + 0.4j, -0.6 + 0.1j, 0.2 + 0.5j, 0.4 - 0.3j]
X1 = [0.0, 1.0 + 0.5j, 0.9 - 0.2j]
X2 = [0.7 - 0.4j, 0.6 + 0.3j]
B1 = [0.3 + 0.2j, 0.5 - 0.1j, 0.4 + 0.6j]
B2 = [0.2 - 0.3j, -0.1 + 0.4j]
B3 = [0.6 + 0.1j, 0.3 - 0.5j]
# Field used by the left drive: the perturbed field that ideal induction would
# give for the same displacement.  It differs from B1..B3 above.
B1_IDEAL = [0.33 + 0.18j, 0.52 - 0.04j, 0.37 + 0.66j]
B2_IDEAL = [0.25 - 0.27j, -0.06 + 0.43j]
B3_IDEAL = [0.55 + 0.16j, 0.36 - 0.47j]
C = 1.7
RNTOR = -2.0
PI2 = 2.0 * math.pi * math.pi


def drive(cell, b1, b2, b3):
    x1 = 0.5 * (X1[cell] + X1[cell + 1])
    q1 = 0.5 * (b1[cell] + b1[cell + 1])
    return H[0] * x1 + H[1] * X2[cell] + H[2] * q1 + H[3] * b2[cell] + H[4] * b3[cell]


def left_coefficient(cell, ideal):
    """KYV2/KYPPARA row that makes the left drive c*conj(h) per cell."""
    b = (B1_IDEAL, B2_IDEAL, B3_IDEAL) if ideal else (B1, B2, B3)
    h = drive(cell, *b)
    return C * h.conjugate() / (X2[cell].conjugate() * G)


def fortran_complex(z):
    return '(%.17e,%.17e)' % (z.real, z.imag)


def array(values):
    return '[' + ','.join(fortran_complex(z) for z in values) + ']'


def driver(kernel, ideal_left):
    d = [left_coefficient(cell, ideal_left) for cell in range(2)]
    return f'''program test_ntv
use globalm
use dimensim
use rcomdm
use kineticm
use torquem
implicit none
complex*16, allocatable :: a(:,:,:),b(:,:,:),c(:,:,:),dm(:,:,:),
& e(:,:,:),f(:,:,:),g(:,:,:),h(:,:,:)
integer js,kgrid,is
NR=2; NRP1=3; MSMAX=1; NCHI=4; NFIT=0; NSPECIES=2; MLMAX=1
NXCOMP=3; NYCOMP=7; MXMAX=NXCOMP*MSMAX; MYMAX=NYCOMP*MSMAX
KXX1=1; KXB1=2; KXV1=3; KYX2=1; KYB2=2; KYB3=3; KYV2=4; KYPR=5
KYPPARA=6; KYPPERP=7
RNTOR={RNTOR}d0; KNTV=21; INCKIN=1; KEFORM=2
ODWKCOM=.true.; ISWEEP=1; NSWEEP=1; ALPHAD=1; PTRAPH=0; CTEDGE=10
INCDPHI=0; IDIAMB=0; IDIAMTE=0; KNTVELL=999
PSPECIES_AP=0; PSPECIES_AT=0; PSPECIES_NP=0; PSPECIES_NTB=0; PSPECIES_NTD=0
PSPECIES_NP(1)=1
allocate(CS(NRP1),CSM(NR),CSH(NR),RM(MSMAX,2),RJAM(NRP1,NCHI))
CS=[.05d0,.4d0,1d0]; CSM=[.2d0,.7d0]; CSH=[.4d0,.6d0]; RM=2; RJAM=1
allocate(X1U(NRP1,1),X2U(NRP1,1),B1U(NRP1,1),B2U(NRP1,1),B3U(NRP1,1))
X1U(:,1)={array(X1)}; B1U(:,1)={array(B1)}
X2U=0; B2U=0; B3U=0
X2U(1:2,1)={array(X2)}; B2U(1:2,1)={array(B2)}; B3U(1:2,1)={array(B3)}
allocate(TORQUENTV(NRP1),TORQUENTVI(NRP1),TORQUENTVE(NRP1))
allocate(HKMIN(NRP1,2),HKMAX(NRP1,2),WFUN(NRP1,2),RLM(1))
HKMIN=0; HKMAX=1; WFUN=1; RLM=1
allocate(VI(4,1,NSPECIES),VPARA(1,1),VPERP(1,1),VDPHI(1,1))
allocate(VX1(1,1),VX2(1,1),VQ1(1,1),VQ2(1,1),VQ3(1,1),VDP(1,1))
VI=0; VI(1,1,1)={fortran_complex(kernel)}
VPARA={fortran_complex(G)}; VPERP=0; VDPHI=0; VDP=0
VX1={fortran_complex(H[0])}; VX2={fortran_complex(H[1])}; VQ1={fortran_complex(H[2])}
VQ2={fortran_complex(H[3])}; VQ3={fortran_complex(H[4])}
allocate(SF0(1,NSPECIES,0:3,4)); SF0=0
call alloc3(VX1PARA); call alloc3(VX1PERP); call alloc3(VX2PARA); call alloc3(VX2PERP)
call alloc3(VQ1PARA); call alloc3(VQ1PERP); call alloc3(VQ2PARA); call alloc3(VQ2PERP)
call alloc3(VQ3PARA); call alloc3(VQ3PERP); call alloc3(VDPPARA); call alloc3(VDPPERP)
call alloc3(VX1PARAM); call alloc3(VX1PERPM); call alloc3(VX2PARAM); call alloc3(VX2PERPM)
call alloc3(VQ1PARAM); call alloc3(VQ1PERPM); call alloc3(VQ2PARAM); call alloc3(VQ2PERPM)
call alloc3(VQ3PARAM); call alloc3(VQ3PERPM); call alloc3(VDPPARAM); call alloc3(VDPPERPM)
call ALLOCATEDWKCOMPMAT
call execute_command_line('mkdir -p DATA_STORAGE')
do kgrid=1,2
do js=1,NRP1
 VX1PARAC=0; VX1PERPC=0; VX2PARAC=0; VX2PERPC=0; VQ1PARAC=0
 VQ1PERPC=0; VQ2PARAC=0; VQ2PERPC=0; VQ3PARAC=0; VQ3PERPC=0
 call KJPFILL(js,js,kgrid,1d0,1,1)
 call WRITE_SURFACE_QUANTITIES(js,kgrid)
enddo
enddo
call DEALLOCATEDWKCOMPMAT
allocate(a(MXMAX,MXMAX,NRP1),b(MXMAX,MXMAX,NRP1),c(MXMAX,MXMAX,NRP1))
allocate(dm(MYMAX,MYMAX,NRP1),e(MXMAX,MYMAX,NRP1),f(MYMAX,MXMAX,NRP1))
allocate(g(MYMAX,MXMAX,NRP1),h(MXMAX,MYMAX,NRP1))
a=0; b=0; c=0; dm=0; e=0; f=0; g=0; h=0
dm(KYV2,KYPPARA,1:2)={array(d)}
! Node-1 radial rows are never contracted (the CALCDWKPROF node row uses I=IS+1).
h(KXV1,KYPPARA,1)=1; h(KXV1,KYPPERP,1)=1
call CALCDWKCOMP(a,b,c,dm,e,f,g,h)
do is=1,NR
 write(*,'(A,I2,ES26.17)') 'TORQUE ',is,TORQUENTVI(is)
enddo
contains
subroutine alloc3(x)
complex*16, allocatable :: x(:,:,:)
allocate(x(MSMAX,MSMAX,NRP1)); x=0
end subroutine
end program
'''


def extract(text, name):
    match = re.search(r'^      SUBROUTINE %s\b' % name, text, re.M)
    if match is None:
        return ''
    end = re.search(r'^      END(\s+SUBROUTINE\s+%s)?\s*$' % name, text[match.start():], re.M)
    return text[match.start():match.start() + end.end()] + '\n'


def run_native(source, build_parent, kernel, ideal_left):
    src = source.resolve() / 'MarsQ_2FK'
    text = (src / 'kinetic.f').read_text()
    with tempfile.TemporaryDirectory(prefix='ntv-torque-', dir=build_parent) as directory:
        work = Path(directory)
        # DERCHI is reached only by the request-file direct check, never here.
        stub = '      SUBROUTINE DERCHI(A,B,N,M)\n      STOP \'UNEXPECTED DERCHI\'\n      END\n'
        (work / 'routines.f').write_text(''.join(extract(text, name) for name in ROUTINES) + stub)
        (work / 'driver.f90').write_text(driver(kernel, ideal_left).replace('\n& ', ' &\n& '))
        names = ['dimensim.f', 'globalm.f', 'rcomdm.f', 'kineticm.f', 'anisotropicm.f',
                 'gijlm.f', 'convolcofm.f', 'torquem.f', 'ToolBox.f90']
        command = ['gfortran', '-O0', '-fdefault-real-8', '-fdefault-double-8',
                   '-ffixed-line-length-none', '-ffree-line-length-none',
                   '-fallow-argument-mismatch', '-ffunction-sections', '-fdata-sections',
                   '-ffpe-trap=zero,invalid',
                   '-fcheck=bounds', '-w', '-I' + str(src)]
        command += [str(src / name) for name in names]
        command += ['routines.f', 'driver.f90', '-Wl,--gc-sections', '-o', 'test']
        compiled = subprocess.run(command, cwd=work, capture_output=True, text=True)
        if compiled.returncode:
            raise RuntimeError(compiled.stdout + compiled.stderr)
        result = subprocess.run([str(work / 'test')], cwd=work, capture_output=True, text=True)
        if result.returncode:
            raise RuntimeError(result.stdout + result.stderr)
        torque = [float(line.split()[2]) for line in result.stdout.splitlines()
                  if line.split()[:1] == ['TORQUE']]
        if len(torque) != 2:
            raise RuntimeError('missing torque rows:\n' + result.stdout[-2000:])
        return torque


def expected_hermitian(kernel):
    # This branch uses T = n Im(dW_k) (upstream issue 21).
    return [RNTOR * PI2 * C * abs(drive(cell, B1, B2, B3)) ** 2 * kernel.imag
            / (4 * math.pi ** 2) for cell in range(2)]


class TestDissipativeNTVTorque(unittest.TestCase):
    source = Path(__file__).resolve().parents[1]
    build_parent = Path.cwd()

    def torque(self, kernel, ideal_left):
        return run_native(self.source, self.build_parent, kernel, ideal_left)

    def test_consistent_drives_give_known_dissipative_torque(self):
        kernel = 0.7 + 0.2j
        for actual, expected in zip(self.torque(kernel, False), expected_hermitian(kernel)):
            self.assertNotEqual(expected, 0.0)
            self.assertAlmostEqual(actual / expected, 1.0, places=12)

    def test_reactive_kernel_gives_zero_torque_with_induction_residual(self):
        scale = max(abs(t) for t in expected_hermitian(0.7j))
        for actual in self.torque(0.7 + 0.0j, True):
            self.assertLessEqual(abs(actual), 1e-14 * scale)

    def test_reactive_part_does_not_change_torque_with_induction_residual(self):
        full = self.torque(0.7 + 0.2j, True)
        dissipative = self.torque(0.2j, True)
        for a, b in zip(full, dissipative):
            self.assertNotEqual(b, 0.0)
            self.assertAlmostEqual(a / b, 1.0, places=12)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--build-parent', type=Path, default=Path.cwd())
    args = parser.parse_args()
    TestDissipativeNTVTorque.source = args.source
    TestDissipativeNTVTorque.build_parent = args.build_parent
    unittest.main(argv=['test'])


if __name__ == '__main__':
    main()
