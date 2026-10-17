"""Compile the executed KIA kernels and compare with a Maxwellian energy oracle."""
import os
from pathlib import Path
import re
import subprocess
import tempfile
import unittest

import numpy as np
from scipy.integrate import quad


ROOT = Path(__file__).resolve().parents[1]


class EnergyUnderflow(unittest.TestCase):
    def test_native_trapped_and_passing_integrals(self):
        source_root = Path(os.environ.get('MARS_ORACLE_SOURCE', ROOT))
        source = (source_root / 'MarsQ_2FK/anisotropic.f').read_text()
        names = ('KIA_TRAP', 'KIA_PASS', 'KENERGYMESH', 'KTRAPDRIFTAT')
        kernels = []
        for name in names:
            start = re.search(rf'^      SUBROUTINE {name}(?:\(|\s*$)', source, re.M).start()
            end = source.index('END SUBROUTINE ' + name, start)
            kernels.append(source[start:end] + 'END SUBROUTINE ' + name + '\n')
        # Derivatives of f(E,psi)=exp(-E)*exp(a psi), with local temperature
        # gradient b: dpsi=f*(a+(E-2.5)*b), dE=-f.  No derivative division.
        def integrand(energy):
            numerator = -2 * energy**2.5 * np.exp(-energy) * (
                -.06 * (.2 + (energy - 2.5) * .3) - (2 + 1j))
            return numerator / (2 + 1j * (1 - .01 * energy**-1.5))
        expected = complex(
            quad(lambda e: integrand(e).real, 0, np.inf, epsabs=1e-11)[0],
            quad(lambda e: integrand(e).imag, 0, np.inf, epsabs=1e-11)[0])
        scratch = ROOT / 'build'
        scratch.mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(dir=scratch, prefix='energy-oracle-') as td:
            work = Path(td)
            (work / 'kernels.f').write_text(''.join(kernels))
            (work / 'oracle.f').write_text('\n'.join('      ' + s for s in DRIVER.splitlines()) + '\n')
            mods = ['globalm.f', 'dimensim.f', 'rcomdm.f', 'kineticm.f', 'anisotropicm.f']
            command = ['gfortran', '-O2', '-fdefault-real-8', '-fdefault-double-8',
                       '-ffixed-line-length-none', '-ffree-line-length-none',
                       '-I' + str(source_root / 'MarsQ_2FK')]
            command += [str(source_root / 'MarsQ_2FK' / f) for f in mods]
            command += ['kernels.f', 'oracle.f', '-o', 'oracle']
            compiled = subprocess.run(command, cwd=work, capture_output=True, text=True)
            self.assertEqual(compiled.returncode, 0, compiled.stdout + compiled.stderr)
            output = subprocess.run(['./oracle'], cwd=work, check=True,
                                    capture_output=True, text=True).stdout
        for line in output.splitlines():
            n, erange, tr, ti, pr, pi = map(float, line.split())
            trapped, passing = complex(tr, ti), complex(pr, pi)
            self.assertTrue(np.isfinite(trapped), f'NEPK{n:g}: nonfinite trapped integral')
            self.assertTrue(np.isfinite(passing), f'NEPK{n:g}: nonfinite passing integral')
            self.assertLess(abs(trapped / expected - 1), 1e-5)
            self.assertLess(abs(passing / expected - 1), 1e-5)
            print(f'NEPK{n:g}: Emax={erange:.9g}; trapped/passing={trapped}; '
                  f'oracle_rel={abs(trapped / expected - 1):.6g}')


DRIVER = '''
program energy_oracle
    use globalm
    use dimensim
    use rcomdm
    use kineticm
    use anisotropicm
    implicit none
    include 'compam.inc'
    integer :: level, k, channel
    complex*16 :: trapped, passing
    integer, parameter :: levels(4) = [101, 201, 401, 801]
    NR=2; NRP1=3; KRADLAYER=0
    allocate(RLM(1), ISPECIES_EK(1), CS(3), T(3), DPSIDS(3))
    allocate(EPSALPHA(3,2,1), NUEFF(3,2,1), OMEGAE0(3,2))
    allocate(DEPSALPHADPSI(3,2,1), SLAM0(1,1))
    RLM=0; ISPECIES_EK=0; CS=[0d0,.5d0,1d0]; T=1; DPSIDS=1
    EPSALPHA=1; NUEFF=.01d0; OMEGAE0=0; DEPSALPHADPSI=0; SLAM0=0
    RNTOR=-3; B0K=1; OMEGACI0=50; OMEGA=dcmplx(-2d0,-1d0)
    NTVNUMEREFAC=1; INUTYPE=1; DRIFT=0; OMEGAB=0; RQK=0
    ISPECIES_F0=0; PSPECIES_FOWP=0; PSPECIES_FOWT=0
    ESPECIES_M=2; ESPECIES_Z=1
    do level=1,4
        NEPK=levels(level); NEPK2=2*NEPK-1
        allocate(EPK(NEPK),ZEPK(NEPK2),ZEPKO(NEPK2),ZEPKN(NEPK2))
        allocate(ZRESU(NEPK2,3,0:6,1))
        call KENERGYMESH
        ZRESU=0
        do channel=1,3
            do k=1,NEPK2
                ZRESU(k,channel,1,1)=exp(-ZEPKN(k))*(.2d0+(ZEPKN(k)-2.5d0)*.3d0)
                ZRESU(k,channel,2,1)=-exp(-ZEPKN(k))
            end do
        end do
        call KIA_TRAP(2,1,1,0,1,.5d0,0,trapped)
        call KIA_PASS(2,1,1,0,1,.5d0,0,passing)
        print *, NEPK, maxval(ZEPKN), real(trapped), aimag(trapped), real(passing), aimag(passing)
        deallocate(EPK,ZEPK,ZEPKO,ZEPKN,ZRESU)
    end do
end program
subroutine PRECOBS_SCALAR
    error stop 'Unexpected precession observer'
end subroutine
subroutine PRECOBS_SAMPLE
    error stop 'Unexpected precession observer'
end subroutine
subroutine PRECOBS_INTEGRAL
    error stop 'Unexpected precession observer'
end subroutine
subroutine KELLTRACESELECT(js,kgrid,otrace)
    implicit none
    integer :: js,kgrid
    logical :: otrace
    otrace=.false.
end subroutine
subroutine WRITEELLTRACEENERGY
    error stop 'Unexpected trace'
end subroutine
subroutine WRITEELLTRACEPITCH
    error stop 'Unexpected trace'
end subroutine
subroutine SPLINE1D
    error stop 'Unexpected pole subtraction'
end subroutine
'''


if __name__ == '__main__':
    unittest.main()
