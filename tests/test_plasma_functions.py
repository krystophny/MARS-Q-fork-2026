"""Analytic energy integrals of the fast kinetic paths (KFASTRUN=1).

KI_PRECESSION evaluates J(W) = int_0^inf x**2.5 exp(-x)/(x+W) dx as
2*sqrt(pi)*(3/8 - W/4 + W**2/2 - W**2*plasmax(W, a)), a = -RTMP.  For W > 0 there is
no resonance and J is real for either sign of a and of a zero imaginary part.
It is compared with scipy quadrature.
"""
import cmath
import math
import re

from scipy.integrate import quad

from support import ROOT, NativeTest

SRC = ROOT / "MarsQ_2FK" / "kinetic.f"


def function(name):
    match = re.search(r"^ +function +%s\b.*?^ +end *$" % name, SRC.read_text(), re.M | re.S | re.I)
    return match.group(0) + "\n"


OMEGAS = (0.5, 2.0, 5.0)
YS = (complex(-3.6, -0.05), complex(-3.6, 0.05), complex(-4.5, -0.1), complex(-2.0, -0.1),
      complex(2.0, 0.3), complex(5.0, 0.5))

DRIVER = """      program plasma_functions
      implicit none
      complex*16 plasmax, plasmay, z, y
      real*8 w(3), a, pi, zero
      complex*16 ys(%d)
      integer i, j, k
      pi = acos(-1.0d0)
      w = [%s]
      ys = [%s]
      do i = 1, 3
         do j = -1, 1, 2
            a = j
            do k = -1, 1, 2
               zero = sign(0.0d0, dble(k))
               z = dcmplx(w(i), zero)
               z = 2*sqrt(pi)*(3./8.-w(i)/4.+w(i)**2/2.-w(i)**2*plasmax(z,a))
               write(*,'(A,3ES24.15,2ES24.15)') 'X', w(i), a, dble(k), z
            end do
         end do
      end do
      do i = 1, size(ys)
         y = plasmay(ys(i))
         write(*,'(A,4ES24.15)') 'Y', ys(i), y
      end do
      end program plasma_functions
"""


class PlasmaFunctions(NativeTest):
    def setUp(self):
        super().setUp()
        main = DRIVER % (len(YS), ", ".join("%.17gd0" % w for w in OMEGAS),
                         ", ".join("(%.17gd0,%.17gd0)" % (y.real, y.imag) for y in YS))
        self.exe = self.compile({"functions.f": function("plasmaz") + function("plasmax") +
                                 function("plasmay"), "main.f90": main})
        self.rows = [line.split() for line in self.execute(self.exe).stdout.splitlines()]

    def test_precession_integral_without_resonance(self):
        rows = [list(map(float, r[1:])) for r in self.rows if r[0] == "X"]
        self.assertEqual(len(rows), 12)
        for w, a, zero, re_, im in rows:
            exact = quad(lambda x: x**2.5*math.exp(-x)/(x+w), 0, math.inf, epsabs=0, epsrel=1e-12)[0]
            with self.subTest(w=w, a=a, zero=zero):
                self.assertLess(abs(complex(re_, im) - exact), 1e-4*exact)


if __name__ == "__main__":
    import unittest
    unittest.main()
