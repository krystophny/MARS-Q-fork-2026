"""BS_ESC response matrix: factorized kernel against the original harmonic-pair loop."""
import os
from pathlib import Path
import re
import subprocess
import tempfile
import unittest

ROOT = Path(os.environ.get("MARS_ORACLE_SOURCE", Path(__file__).resolve().parents[1]))

ORACLE = """program esc_response_oracle
implicit none
integer, parameter :: nchi = 48, np = 24, nm = 21, mlo = 4, mhi = 18
real*8 :: rmode(nm), rntor, hchi, pi, t, temp1, temp2
real*8 :: r1(nchi), z1(nchi), r2(nchi), z2(nchi), drdc(nchi), dzdc(nchi), dsdr(nchi), dsdz(nchi)
real*8 :: pc(np), ps(np), pcn(np), psn(np), hphi
complex*16 :: eca(nm, nm), ref(nm, nm), ctmp1, ctmp2
complex*16, parameter :: ci = (0d0, 1d0)
integer :: j, j1, j2, ms, ks
pi = 4d0*atan(1d0); hchi = 2*pi/nchi; hphi = 2*pi/np; rntor = -3d0
do j = 1, nm
  rmode(j) = j - 11
end do
! IFEED (1) outside IBNM (2), both shaped and shifted, as in BS_ESC
do j = 1, nchi
  t = hchi*(j - 1)
  r1(j) = 1.0d0 + 0.62d0*cos(t + 0.21d0*sin(t)); z1(j) = 1.15d0*sin(t)
  r2(j) = 1.0d0 + 0.41d0*cos(t + 0.13d0*sin(t)); z2(j) = 0.03d0 + 0.78d0*sin(t)
  drdc(j) = -0.62d0*sin(t + 0.21d0*sin(t))*(1 + 0.21d0*cos(t)); dzdc(j) = 1.15d0*cos(t)
  dsdr(j) = 0.78d0*cos(t)*r2(j)
  dsdz(j) = 0.41d0*sin(t + 0.13d0*sin(t))*(1 + 0.13d0*cos(t))*r2(j)
end do
do j = 1, np
  temp1 = hphi*(j - 1) - pi
  pc(j) = cos(temp1); ps(j) = sin(temp1); pcn(j) = cos(rntor*temp1); psn(j) = sin(rntor*temp1)
end do
! original BS_ESC loop (marsq.f before the factorization)
ref = (0d0, 0d0)
do ms = mlo, mhi
do ks = mlo, mhi
do j1 = 1, nchi
do j2 = 1, nchi
  ctmp1 = exp((rmode(ks)*(j1-1) - rmode(ms)*(j2-1))*hchi*ci)
  temp1 = (r2(j2)-r1(j1))**2 + (z2(j2)-z1(j1))**2
  do j = 1, np
    temp2 = temp1 + 2.*r1(j1)*r2(j2)*(1.-pc(j))
    temp2 = temp2*sqrt(temp2)
    ctmp2 = (dsdr(j2)*(r1(j1)*(z1(j1)-z2(j2))*pc(j)*pcn(j)*rmode(ks)/rntor + &
             (r1(j1)*dzdc(j1)+(z2(j2)-z1(j1))*drdc(j1))*ps(j)*psn(j)*ci) + &
            dsdz(j2)*(r1(j1)*pcn(j)*(-r1(j1)+r2(j2)*pc(j))*rmode(ks)/rntor - &
             r2(j2)*ps(j)*psn(j)*drdc(j1)*ci))*ctmp1/temp2
    ref(ms, ks) = ref(ms, ks) + ctmp2
  end do
end do
end do
end do
end do
eca = (7d0, 7d0)
call esc_response(nchi, np, nm, mlo, mhi, rmode, rntor, r1, z1, r2, z2, drdc, dzdc, dsdr, dsdz, &
                  pc, ps, pcn, psn, hchi, eca)
if (maxval(abs(eca - ref)) > 1d-12*maxval(abs(ref))) then
  print *, maxval(abs(eca - ref)), maxval(abs(ref))
  error stop 'factorized ESC response differs from the harmonic-pair loop'
end if
if (any(eca(1:mlo-1, :) /= 0) .or. any(eca(:, mhi+1:) /= 0)) error stop 'entries outside the ESC range set'
print *, 'max relative difference', maxval(abs(eca - ref))/maxval(abs(ref))
end program
"""


class EscResponse(unittest.TestCase):
    def test_factorized_kernel_matches_original_loop(self):
        text = (ROOT / "MarsQ_2FK" / "marsq.f").read_text()
        m = re.search(r"^      SUBROUTINE ESC_RESPONSE\(.*?^      END *$", text, re.M | re.S)
        self.assertIsNotNone(m, "ESC_RESPONSE not found in marsq.f")
        with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as tmp:
            work = Path(tmp)
            (work / "esc.f").write_text(m.group(0) + "\n")
            (work / "oracle.f90").write_text(ORACLE)
            fc = os.environ.get("FC", "gfortran")
            subprocess.run([fc, "-O0", "-fcheck=all", "-ffixed-line-length-none", "esc.f", "oracle.f90", "-o", "oracle"],
                           cwd=work, check=True, capture_output=True, text=True)
            run = subprocess.run(["./oracle"], cwd=work, capture_output=True, text=True)
            self.assertEqual(run.returncode, 0, run.stdout + run.stderr)


if __name__ == "__main__":
    unittest.main()
