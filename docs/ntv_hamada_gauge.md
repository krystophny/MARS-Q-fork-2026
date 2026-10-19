# |B| in perturbed Hamada coordinates for the NTV module

`torque.f` calls `ntv_predata` (`NTV_pre::MARS_out`) and `gntv1d`, the
Shaing-type NTV torque in Hamada coordinates `(V, theta_H, zeta_H)`. Its input
`bmn` is the spectrum of the field-strength perturbation. This note derives
that input and the implementation in `cal_xi_liu`.

Evidence labels: **[proved]** means a FortSym native zero test returned true.
The script is [`tools/ntv_hamada_gauge/identities.py`](../tools/ntv_hamada_gauge/identities.py)
and its verdicts are in [`fortsym_verdicts.json`](../tools/ntv_hamada_gauge/fortsym_verdicts.json).
**[E\*, proved]** cites an identity proved in itpplasma/rmp_torque c3614334,
`docs/boozer_gauge_ntv_derivation.md` (PR #55). **[numeric]** names a test in
`tests/test_ntv_hamada_gauge.py`. **[cited]** marks a paper equation, and
**[hypothesis]** marks a statement that is not verified.

## 1. What was wrong

Before this change, `cal_xi_liu` built

    dB_old = xi.grad B,   xi = xi_perp + X3 B,

where `X3 = xi.B/B^2` is MARS-F's fluid parallel displacement
(`ntv_pre.f` `cal_metric`, `gg3 = B^chi`). This has two defects.

1. `X3` is the MHD parallel velocity response. It is not the parallel
   displacement of the Hamada coordinate gauge. A purely field-aligned `X3 B`
   leaves the field unchanged (`curl(lam B x B) = 0`, E2.1x), but it changed
   `dB_old` by `X3 B.grad B`.
2. The Eulerian part `b.deltaB` of the Lagrangian field strength is missing.
   It is `dB_L = b.deltaB + xi.grad B` (Park 2009 Eq. 7 **[cited]**).

## 2. Perturbed Hamada coordinates

The displaced coordinates have the Jacobian `J' = J(1 + div xi)`
**[E4.1, proved]**. A Hamada base has `J` constant on each surface, so the
perturbed coordinates stay Hamada iff `div xi = c(psi)`. For a single
`n != 0`, `c = 0`. With `xi = xi_perp + lam B` and `div B = 0`,
`div(lam B) = B.grad lam`. This gives the magnetic differential equation

    (2.1)  B.grad lam = c - div xi_perp,     c = 0 for n != 0.

This matches Park 2009 Eq. 10 **[cited]** and rmp_torque #55 Section 4
("Hamada and divergence-free gauges", (4.13)). The sign has been checked
against both. `lam = xi_par/B` is fixed up to a flux function. That flux
function vanishes for `n != 0` off resonance.

The field strength at fixed perturbed Hamada labels is the Lagrangian
`dB_L(xi)`. Using E1.3 and E1.5 **[proved]**:

    (2.2)  dB_H = dB_L(xi_perp) + lam B.grad B = -B(d + k) + lam B.grad B,
           d = div xi_perp,  k = kappa.xi_perp.

Equation (2.2) holds under ideal induction, `deltaB = curl(xi x B)`, which is
the field-line-preserving map that defines the perturbed coordinates. In
force balance,

    (2.3)  k = xi_perp.grad B/B + mu0 p'(psi) psi'(s) X1/B^2              [H5 proved]

This is the form MARS-K uses (`RX1P = 2 PPEQ J/B`). It needs `PPEQ`
(`mu0 dp/dpsi`), not second derivatives of `R, Z`.

## 3. MARS chart

MARS uses `(s, chi, phi)` with equal-arc `chi` and `Jac = R(R_s Z_chi - R_chi Z_s)`
(`cal_metric`). It also uses `B^chi = psi'/J`, `B^phi = F/R^2` and `B_phi = F`.
The displacement is (`cal_metric` comment, controller-verified)

    xi = X1 (e_s - B_s B/B^2) + X2 J (B x grad s)/B^2 + X3 B,

with contravariant components

    xi^s = X1,  xi^chi = X1 gg1 + X2 gg2,  xi^phi = -X1 B_s B^phi/B^2 - X2 B_chi/B^2,
    xi_perp.B = 0                                                                 [H1 proved]

Then, with `D_s` and `D_chi` acting on the MARS harmonics and `d_phi -> i n`:

    J d = dm1 X1 + J D_s X1 + dm2 D_chi X1 + dm3 X2 + dm4 D_chi X2
          + i n (dm5 X1 + dm6 X2)                                                  [H2 proved]
    dm2 = -g_schi psi'^2/(J B^2),  dm4 = J F/B^2,
    dm5 = -g_schi psi' F/(R B)^2,  dm6 = -g_chichi psi'/B^2,
    dm1 = d_s J + d_chi dm2,       dm3 = d_chi dm4.

`cal_divmetric` computes these coefficients on the `(s, chi)` grid, using
periodic centred differences in `chi`. They are mapped to Hamada angles with
the other metric arrays.

**Map to Hamada angles.** `coord_trans` defines
`theta_H = 2 pi int J dchi/V_S` and `zeta_H = phi + dphi(chi)`, with
`dphi' = q(1 - V_S/(R^2 K)) theta_H'`. Here `V_S = oint J dchi` and
`K = oint J/R^2 dchi`. Then

    B.grad theta_H = 2 pi psi'/V_S = B^theta (a flux function),                    [H3a proved]
    B.grad zeta_H  = q B^theta,   provided F = 2 pi q psi'/K.                      [H3b proved]

The proviso is that the profile `q` is the straight-field-line `q` of the
chart. `MARS_out` prints `max |F K/(2 pi psi' q) - 1|` on every call, so that
this can be checked on a production equilibrium **[hypothesis until a
qualification run prints it]**. In this chart `B.grad = B^theta(d_theta + q d_zeta)`.
`J_H` is constant, so `(theta_H, zeta_H)` are Hamada angles.

`trans_1d` evaluates `sum_m X_m exp(i(m chi - n dphi))`, which is the
perturbation `exp(i n phi)` on `zeta_H = 0`. The `bmn` projection on the
`theta_H` midpoint grid therefore returns the coefficients of
`exp(i(k theta_H + n zeta_H))`. `gntv1d` receives these coefficients for the
single `n = RNTOR`, as before.

## 4. Fourier solution and regularisation

On `exp(i(k theta_H + n zeta_H))`, equation (2.1) gives

    (4.1)  lam_k = i d_k/(B^theta (k + nq))                                          [H4a proved]

with the Park 2009 regularisation **[cited]**

    (4.2)  lam_k = i d_k (k + nq)/(B^theta ((k + nq)^2 + sigma^2)).

The residual of (4.2) is `d_k sigma^2/((k+nq)^2 + sigma^2)` **[H4b proved]**.
It is non-zero only inside the layer `|k + nq| <~ sigma`. The code uses one
radial cell, `sigma = |n dq/ds| ds`: centred `(s_{i+1} - s_{i-1})/2`, and
`s_i - s_{i-1}` at the edge. It sets `lam_k = 0` when `(k + nq)^2 + sigma^2 = 0`;
this case is `n = 0, k = 0` and is the constant `c`. The cell width is a choice
**[hypothesis]**: Park reports the torque insensitive to `sigma` within a factor
of 10, and #55 Section 6 gives the physical layer width. The DFT of `d` uses
all `n_chi` harmonics of the `theta_H` grid. `lam` is synthesised on the same
grid.

Implemented spectrum (`cal_xi_liu`):

    dB_H = -B d - xi_perp.grad B - mu0 p' psi' X1/B + lam B^chi d_chi B.

`X3` is read by `get_disp` but not used.

## 5. Verification

- FortSym (`identities.py`): H1-H5 above all **[proved]**, as are the MARS-K
  identities H6-H8 (Section 6).
- **[numeric]** (a) A field-aligned `X3` alone gives `bmn = 0` exactly. Adding
  `X3` to any `X1, X2` leaves `bmn` bitwise unchanged.
- **[numeric]** (b) On a Hamada base with `B0 = b00(1 - eps cos theta)`,
  `div xi_perp != 0`, `cp != 0` and `sigma = 0.2`, `bmn` equals rmp_torque's
  `delta_b_boozer(dB_L, source = div xi_perp)` to `1e-11`. That function is the
  #55 kernel; with source `T = div xi_perp`, `w = lam/J_H` solves `D w = -d`,
  which is the Hamada gauge (Section 2). The reference is generated by
  `tools/ntv_hamada_gauge/rmp_torque_reference.py`.
- **[numeric]** (c) The circular case uses `eps = 0.1`, `q = 1.1 + 0.8 s^2`,
  `n = 2` and five m-harmonics of `X1` and `X2`, and runs the full chain
  `cal_metric` -> `cal_divmetric` -> `coord_trans` -> `coord_map` -> `cal_xi_liu`.
  The independent oracle uses an exact Hamada chart
  (`theta_H = chi + e sin chi`, closed-form `dphi`) and Cartesian finite
  differences of the explicit vector `xi_perp`, of `|B|` and of `B.grad|B|`.
  It also checks that `B.grad lam = -d` along the Cartesian field line, to
  `1e-5`. This check fails if the sign of `k + nq` is flipped. The relative
  agreement of the spectra is better than `1e-3` on three surfaces.

## 6. MARS-K kinetic drive

`kinetic.f` (around lines 9017-9026 at 12c9809/55739a9; coefficients from line
6900) implements Liu et al 2008 Eq. 11 **[cited]**:
`H_L = m v_par^2 kappa.xi_perp + mu(Q_L,par + grad B.xi_perp)`. Normalised by
`E` with `z = mu B/E`, this is the ledger's

    h_M = z b_E + (2 - z) a_perp + 2(1 - z) c_p.

- With force balance (H5), `h_M = z(b_E + a_perp) + 2(1 - z) k = h_id`
  **[H6a proved]**. With ideal induction, `h_id = -z d + (2 - 3z) k`
  **[H6b proved]**.
- The Logan 2013 Eq. 5 / Park 2011 Eqs. 10-11 weight is
  `-[(2E - 3 mu B) dB_L(xi)/B + (2E - 2 mu B) div xi]`. On `xi_perp` it equals
  `E h_id` **[H8 proved]**. On a field-aligned `lam B` it equals
  `-v d_l(m v lam B)`, a total streaming derivative **[H7 proved]**. The
  Jacobian term `(2E - 2 mu B) div xi` is what removes the `xi_par` dependence.
  MARS-K never uses `xi_par` (`X3` does not enter `RVAK*`), so it needs no
  Jacobian term.
- The `v_par^2` and `v^2` split: the `(1 - z)` weight on `c_p` is
  `m v_par^2/(2E)` and multiplies only the pressure part of the curvature. The
  `(2 - z)` weight on `a_perp` is `(m v_par^2 + mu B)/E`. Both follow from H6a.
- The ledger residual is `h_N - h_M = 2(1 - z)(b_E - c_p)` **[H6c proved]**.
  It equals `2(1 - z) S` with `S = -(d + 2k)` **[H6d proved]**. This is #55's
  E2.8b. It is the error of NEO-RT driven by the `xi_perp`-only `dB_L`, not a
  MARS-K defect.

Verdict: MARS-K's drive is gauge-free and equals the gauge-invariant `h_id`
under force balance and ideal induction. It needs no correction. In a
resistive layer, `Q` differs from `curl(xi x B)`. There MARS-K uses the
actual field, which is the physical choice, and the H6b form does not apply.

## References

- Park, Boozer, Menard, Phys. Plasmas 16, 056115 (2009), Eqs. 7, 10.
- Shaing, Phys. Plasmas 10, 1443 (2003), Eqs. 3-5, A3.
- Boozer and Nuehrenberg (2006), Eq. 42.
- Liu et al, Phys. Plasmas 15, 112503 (2008), Eq. 11.
- Park (2011), Eqs. 10-11; Logan et al (2013), Eq. 5.
