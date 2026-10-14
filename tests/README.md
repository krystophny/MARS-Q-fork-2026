# MARS regression tests

The fast standard-library test suite covers the complete local delta from the
upstream baseline `8824bb1`:

| Change | Fast coverage | System coverage |
|---|---|---|
| portable namelist errors | source contract and black-box malformed input | compiler matrix |
| strict `INPUT_BNM` header | parser contract | MAST-U `KKF=-3` run |
| OpenMP equivalent-current kernel | scoping and arithmetic contract | 1-thread/parallel `CURHARMO` comparison |
| GNU/ifx/NVHPC targets | target/flag contract | clean build matrix |
| NVHPC `-O1` workaround | exact flag contract | NVHPC build |
| CHEASE legacy initialization | GNU/ifx flag contract | ITER EQDSK-to-OUTRMAR run |
| perturbative MARS-K NTV | input contracts and unsupported cross-term rejection | paired MARS-F/MARS-K field comparison |
| external frozen B/X import | strict/atomic parser contract and invalid-mode rejection | exact BPLASMA/XPLASMA round trip plus MARS-K torque |
| perturbative passive-operator ordering | every final-sweep `KNTV=21` field, including ordinary native `IPERTURB=0`, assembles passive `KPBKEY=1` blocks once before torque; only imported B/X is restored | paired native/reload B/X and DWK equality plus matching nonzero `TORQUENTV.OUT` |
| MARS-K matrix call contract | all `CALCDWKCOMP` callers pass eight live matrices | debug cached continuation reaches all DWK and torque outputs |
| validated DWK cache recovery | narrow namelist contract and KJP bypass | failed-after-cache ITER runs resume without coefficient regeneration |
| DWK Maxwellian tail | compiled native passing/trapped kernels, finite-birth partners, independent analytic Maxwellian quadrature at NEPK=201/401/801, checked bounds | NEPK=401 frozen-field smoke and NEPK=201 torque regression |
| deterministic frequency diagnostic | `RUU2` is thread-private with its per-surface consumers | repeated OpenMP `FREQUENCIES.OUT` comparison |
| two-species frequency diagnostic | hot-ion `SLAM0(:,3)` is read only when species three exists | GNU `-fcheck=all` MAST-U run |
| retained DWK component workspace | idempotent allocation, component-map, and shape contracts | KNTV=21 reaches `TORQUENTV.OUT` after kinetic assembly |
| independent DWK density check | numerical angular-quadrature oracle and source normalization contract | opt-in direct check closes against the component density without an extra `4*pi^2` factor |
| executable build provenance | profile/flag and manifest hash/tamper contracts | clean-tree build plus independent manifest verification |

Run fast tests with:

```sh
make -C MarsQ_2FK test
```

Add executable-level tests with:

```sh
make -C MarsQ_2FK test-runtime MARS_EXE="$PWD/build/marsq-gnu.x"
```

The long MAST-U and ITER system runs use external/private fixtures and are
recorded in the consuming project rather than copied into this public source
repository.
