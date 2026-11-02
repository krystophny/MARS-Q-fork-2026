"""Temperature in the Shaing torque prefactor of each species (NTV_pre::cal_nu).

gntv scales the torque density with T0 ~ rho_i*q**2*R0**2*wti**2, taking
rho_i = nu_prof(:,4) and wti = nu_prof(:,6) from cal_nu.  For the torque of
a species this is 2 n T of that species, so at equal density the electron
prefactor (Z = -1) divided by the ion prefactor (Z = 1) must be T_e/T_i.
"""
import re

from support import ROOT, NativeTest

SRC = ROOT / "MarsQ_2FK" / "ntv_pre.f"
LIB = ROOT / "MarsQ_2FK" / "ntv_lib.f"


def routine(name, source=SRC):
    match = re.search(r"^ +subroutine +%s\b.*?^ +end subroutine +%s\b" % (name, name),
                      source.read_text(), re.M | re.S | re.I)
    if match is None:
        raise ValueError("Native routine not found: " + name)
    return match.group(0) + "\n"


PROGRAM = """      program electron_temperature
      implicit none
      integer, parameter :: nr = 4, np = 15
      real :: r(nr), drdpsi(nr), rp(nr,np), oi(nr,np), oe(nr,np)
      integer :: i
      r = [0.2, 0.4, 0.6, 0.8]
      drdpsi = [0.9, 0.7, 0.5, 0.4]
      rp = 0.0
      rp(:,1) = r
      rp(:,2) = [1.1, 1.9, 3.2, 4.7]
      rp(:,3) = [0.9, 0.8, 0.6, 0.4]
      rp(:,4) = [9.0, 7.0, 4.0, 1.5]
      rp(:,5) = [6.0, 6.5, 2.0, 1.2]
      rp(:,7) = [-0.1, -0.3, -0.6, -1.1]
      rp(:,8) = [-3.0, -6.0, -9.0, -7.0]
      rp(:,9) = [-2.5, -5.0, -8.0, -6.0]
      rp(:,10) = 0.3*r
      rp(:,13) = -1.0
      rp(:,14) = -1.0
      call cal_nu(0, nr, r, drdpsi, np, rp, 2.0, 6.2, 5.3, 1, 1, 1.17, oi)
      call cal_nu(0, nr, r, drdpsi, np, rp, 2.0, 6.2, 5.3, -1, 1, 1.17, oe)
      do i = 1, nr
         write(*,'(3es24.16)') rp(i,4)/rp(i,5), oi(i,4)*oi(i,6)**2,
     &        oe(i,4)*oe(i,6)**2
      end do
      end program electron_temperature
"""


class NtvElectronTemperature(NativeTest):
    def test_electron_prefactor_scales_with_electron_temperature(self):
        exe = self.compile({"cal_nu.f": routine("cal_nu"),
                            "cald_1d.f": routine("cald_1d", LIB), "main.f": PROGRAM})
        rows = [list(map(float, line.split()))
                for line in self.execute(exe).stdout.split("\n") if line.strip()]
        self.assertEqual(len(rows), 4)
        for te_ti, ion, electron in rows:
            self.assertGreater(ion, 0.0)
            self.assertAlmostEqual(electron / ion, te_ti, delta=1.0e-12 * te_ti)


if __name__ == "__main__":
    import unittest
    unittest.main()
