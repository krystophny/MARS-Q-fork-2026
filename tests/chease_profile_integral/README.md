# Native CHEASE cubic-profile integral gate

- Scope: `CheaseMerge/chease.f::ISOFUN`, both `NTMF0` integration branches. `TMF` holds F²/2 during integration and becomes sqrt(2*TMF) afterwards. Native `SPLINE` defines `D2TTP=d²(TTP)/dψ²`.
- Exact cell integral: h*(0.5*(f0+f1) − h²*(M0+M1)/24). Endpoint direction changes the accumulation sign, not this integral.
- Baseline uses h³/48 in both branches and adds the curvature correction in the upward branch. Constants/linear profiles are unaffected; curved source laws change reconstructed F.
- Extract literal selected production branches and native SPLINE/TRIDAG into the build directory. Original source remains the authority; no hand-written replacement native routine or CAS dependency enters upstream.
- Oracle: exact polynomial antiderivatives, degrees0–3 on seven nonuniform negative-to-zero flux knots. Independently verify native spline second derivatives. Test both integration directions and reject four lost-curvature mutations.
- Baseline8/12; repaired12/12. Four quadratic/cubic failures disappear; constants/linear results remain unchanged. This certifies the selected profile primitive, not full equilibrium/force or E0 edge-q.
- Reproduce: `cmake -S tests/chease_profile_integral -B /tmp/chease-profile-gate`; `cmake --build /tmp/chease-profile-gate`; `ctest --test-dir /tmp/chease-profile-gate --output-on-failure`.
- GNU native arithmetic flags preserve default real8/double8. Full legacy build additionally needs `-mcmodel=medium`; the small model overflows existing large COMMON relocations. Existing accepted binaries/runs remain unchanged.
