"""FortSym zero tests for docs/ntv_hamada_gauge.md.

Usage (FortSym checkout with its C ABI library built)::

    PYTHONPATH=$FORTSYM/python python3 tools/ntv_hamada_gauge/identities.py

Writes ``fortsym_verdicts.json`` next to this file.  A check is PROVED only
when FortSym's native zero decision returns True; probes are not accepted.
Labels H* refer to the derivation note.  Labels E* cited there are proved in
itpplasma/rmp_torque c3614334 (docs/boozer_gauge_ntv_derivation.md).
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import fortsym
import fortsym.sympy as sp

HERE = Path(__file__).resolve().parent


def mars_chart():
    """H1-H3: displacement components, div xi_perp expansion, Hamada chart."""
    s, c, p = sp.symbols("s chi phi")
    J, R, gsc, gcc, B = (sp.Function(n)(s, c) for n in ("J", "R", "g_schi", "g_chichi", "B"))
    psip, F = sp.Function("psip")(s), sp.Function("F")(s)
    X1, X2 = sp.Function("X1")(s, c, p), sp.Function("X2")(s, c, p)
    Bc, Bp = psip / J, F / R**2                     # B^chi, B^phi
    Bs_, Bc_, Bp_ = gsc * Bc, gcc * Bc, F           # covariant B_s, B_chi, B_phi
    B2 = Bc * Bc_ + Bp * Bp_
    # xi_perp = X1 (e_s - B_s B/B^2) + X2 J (B x grad s)/B^2; Levi-Civita gives
    # (J B x grad s)^i = (0, B_phi, -B_chi).  Here |B|^2 is the metric B2.
    xi = [X1, X1 * (-Bs_ * Bc / B2) + X2 * Bp_ / B2, X1 * (-Bs_ * Bp / B2) - X2 * Bc_ / B2]
    gg1, gg2 = -gsc * Bc**2 / B2, F / B2            # cal_metric gg(:,:,1:2)
    out = [("H1a", "xi_perp.B = 0", xi[1] * Bc_ + xi[2] * Bp_ + xi[0] * Bs_),
           ("H1b", "xi^chi = X1 gg1 + X2 gg2", xi[1] - X1 * gg1 - X2 * gg2)]
    # H2 with an arbitrary function B in place of |B| (the expansion is algebraic).
    xiB = [X1, X1 * (-gsc * Bc**2 / B**2) + X2 * F / B**2,
           X1 * (-gsc * Bc * Bp / B**2) - X2 * gcc * Bc / B**2]
    dm2 = -gsc * psip**2 / (J * B**2)
    dm4 = J * F / B**2
    dm5 = -gsc * psip * F / (R * B)**2
    dm6 = -gcc * psip / B**2
    dm1 = J.diff(s) + dm2.diff(c)
    dm3 = dm4.diff(c)
    out.append(("H2a", "J xi^chi = dm2 X1 + dm4 X2", J * xiB[1] - dm2 * X1 - dm4 * X2))
    out.append(("H2b", "J xi^phi = dm5 X1 + dm6 X2", J * xiB[2] - dm5 * X1 - dm6 * X2))
    div = ((J * xiB[0]).diff(s) + (J * xiB[1]).diff(c) + (J * xiB[2]).diff(p)) / J
    code = (dm1 * X1 + J * X1.diff(s) + dm2 * X1.diff(c) + dm3 * X2 + dm4 * X2.diff(c)
            + dm5 * X1.diff(p) + dm6 * X2.diff(p)) / J
    out.append(("H2c", "div xi_perp = cal_xi_liu divx expansion", div - code))
    # H3: coord_trans chart.  VS = oint J dchi, K = oint J/R^2 dchi; q is
    # straight-field-line iff F = 2 pi q psi'/K (the printed chart check).
    Jp, Rp, ps, q, VS, K = sp.symbols("J R psip q V_S K", positive=True)
    Fq = 2 * sp.pi * q * ps / K
    dth = 2 * sp.pi * Jp / VS                       # d theta_H/d chi (theHc)
    ddphi = q * (1 - VS / (Rp**2 * K)) * dth        # d dphi/d chi (fac)
    bth = 2 * sp.pi * ps / VS
    out.append(("H3a", "B.grad theta_H = 2 pi psi'/V_S = bth", (ps / Jp) * dth - bth))
    out.append(("H3b", "B.grad zeta_H = q B^theta",
                Fq / Rp**2 + (ps / Jp) * ddphi - q * bth))
    return out


def gauge_solution():
    """H4: Fourier solution of B.grad lam = -d and its regularisation."""
    k, n, q, bth, sig = sp.symbols("k n q B_theta sigma", positive=True)
    d = sp.Symbol("d_k")
    x = k + n * q
    lam = sp.I * d / (bth * x)
    lam_s = sp.I * d * x / (bth * (x**2 + sig**2))
    return [("H4a", "i B^theta (k+nq) lam_k = -d_k", sp.I * bth * x * lam + d),
            ("H4b", "regularised residual = d_k sigma^2/((k+nq)^2+sigma^2)",
             sp.I * bth * x * lam_s + d - d * sig**2 / (x**2 + sig**2))]


def _jet():
    x, y, z = sp.symbols("x y z")
    X = (x, y, z)
    beta = sp.Symbol("beta", positive=True)
    B = [[0, 0, beta][i] + sum(sp.Symbol(f"g_{i}{j}") * X[j] for j in range(3))
         for i in range(3)]
    a = [sp.Symbol(f"a_{i}") + sum(sp.Symbol(f"a_{i}{j}") * X[j] for j in range(3))
         for i in range(3)]
    return X, B, a


def curvature_force_balance():
    """H5: kappa.xi_perp = xi_perp.grad B/B + xi_perp.[(curl B) x B]/B^2 on 1-jets."""
    X, B, a = _jet()
    d = lambda f, c: f.diff(c)
    cross = lambda u, v: [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                          u[0] * v[1] - u[1] * v[0]]
    dot = lambda u, v: sum(u[i] * v[i] for i in range(3))
    Bm = sp.sqrt(dot(B, B))
    b = [c / Bm for c in B]
    kappa = [sum(b[j] * d(b[i], X[j]) for j in range(3)) for i in range(3)]
    curl = [d(B[2], X[1]) - d(B[1], X[2]), d(B[0], X[2]) - d(B[2], X[0]),
            d(B[1], X[0]) - d(B[0], X[1])]
    xi = cross(B, a)
    res = (dot(kappa, xi) - sum(xi[j] * d(Bm, X[j]) for j in range(3)) / Bm
           - dot(xi, cross(curl, B)) / Bm**2)
    for c in X:
        res = res.subs(c, 0)
    return [("H5", "kappa.xi_perp = xi_perp.grad B/B + xi_perp.(curl B x B)/B^2", res)]


def kinetic_drive():
    """H6-H8: MARS-K source, gauge-invariant action and the Logan/Park weight."""
    z, bE, a, c, d = sp.symbols("z b_E a_perp c_p d")
    k = a + c                                       # H5 with force balance
    bE_ideal = -(d + k) - a                         # E1.2/E1.3 (ideal induction)
    h_M = z * bE + (2 - z) * a + 2 * (1 - z) * c    # Liu 2008 Eq. 11 / kinetic.f
    h_id = z * (bE + a) + 2 * (1 - z) * k
    h_N = (2 - z) * (bE + a)
    S = -(d + 2 * k)
    out = [("H6a", "h_M = z(b_E+a) + 2(1-z)kappa.xi (force balance)", h_M - h_id),
           ("H6b", "h_id = -z d + (2-3z)k (ideal induction)",
            h_id.subs(bE, bE_ideal) - (-z * d + (2 - 3 * z) * k)),
           ("H6c", "h_N - h_M = 2(1-z)(b_E - c_p)", h_N - h_M - 2 * (1 - z) * (bE - c)),
           ("H6d", "b_E - c_p = S = -(d + 2k) (ideal)", (bE - c - S).subs(bE, bE_ideal))]
    l = sp.Symbol("l")
    E, mu, m = sp.symbols("E mu m", positive=True)
    Bl, lam = sp.Function("B")(l), sp.Function("lam")(l)
    v = sp.sqrt(2 * (E - mu * Bl) / m)
    # Logan 2013 Eq. 5 weight -[(2E-3muB) dB_L(xi)/B + (2E-2muB) div xi], field-
    # aligned part xi = lam B: dB_L = lam B B', div xi = B lam' (B.grad = B d_l).
    W_par = -((2 * E - 3 * mu * Bl) * lam * Bl.diff(l) + (2 * E - 2 * mu * Bl) * Bl * lam.diff(l))
    out.append(("H7", "Logan weight of lam B = -v d_l(m v lam B) (total d/dt)",
                W_par + v * (m * v * lam * Bl).diff(l)))
    dd, kk = sp.symbols("d k")
    Bs = sp.Symbol("B", positive=True)
    W_perp = -((2 * E - 3 * mu * Bs) * (-(dd + kk)) + (2 * E - 2 * mu * Bs) * dd)
    zz = mu * Bs / E
    out.append(("H8", "Logan weight of xi_perp = E h_id",
                W_perp - E * (-zz * dd + (2 - 3 * zz) * kk)))
    return out


def all_identities():
    return mars_chart() + gauge_solution() + curvature_force_balance() + kinetic_drive()


def main() -> int:
    rows, failed = [], []
    for label, equation, residual in all_identities():
        value = residual.is_zero
        verdict = {True: "PROVED", False: "NONZERO"}.get(value, "UNKNOWN")
        rows.append({"label": label, "equation": equation, "verdict": verdict})
        if verdict != "PROVED":
            failed.append(label)
        print(f"{label:5s} {verdict:8s} {equation}")
    root = Path(fortsym.__file__).resolve().parents[2]
    rev = subprocess.run(["git", "-C", str(root), "rev-parse", "HEAD"],
                         capture_output=True, text=True).stdout.strip() or "unknown"
    (HERE / "fortsym_verdicts.json").write_text(json.dumps(
        {"engine": "fortsym native zero test", "fortsym_revision": rev,
         "checks": rows, "failed": failed}, indent=2) + "\n")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
