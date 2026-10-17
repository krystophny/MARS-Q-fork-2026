"""Native maintained-background profiles and derivative/normalization oracles."""
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
N = 9
MASS_UNIT = 1.67261e-27
MASSES = (2 * MASS_UNIT, MASS_UNIT / 1836.15267343, 3 * MASS_UNIT)
DENSITIES = (1e19, 2e19, 1e19)
ENERGIES = (1.3e-16, 3.7e-16, 2.1e-16)
NREF = 5e19
PREF = 2.3**2 / (4e-7 * math.pi)


def routine(source, name):
    match = re.search(r"^      subroutine +" + name + r"(?:\(|\s*$)", source,
                      re.M | re.I)
    if match is None:
        raise ValueError("Native subroutine missing: " + name)
    end = re.search(r"^      end(?: subroutine +" + name + r")?\s*$",
                    source[match.start():], re.M | re.I)
    return source[match.start():match.start() + end.end()] + "\n"


def profile(slope=0, temperature_slope=0):
    records = ["MARS_KINETIC_PROFILES_V1",
               "MAINTAINED_FROZEN_FIELD NATIVE_CS DENSITY_M3 TEMPERATURE_ENERGY_J",
               f"3 {N}"]
    records += [f"{k + 1} {z} {mass:.17e}"
                for k, (z, mass) in enumerate(zip((1, -1, 1), MASSES))]
    for j in range(N):
        psi = j / (N - 1)
        values = [math.sqrt(psi)]
        for density, energy in zip(DENSITIES, ENERGIES):
            values += [density * (1 + slope * psi), energy * (1 + temperature_slope * psi**2)]
        records.append(" ".join(f"{value:.17e}" for value in values))
    return records


class ExplicitKineticProfiles(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = Path(os.environ.get("MARS_ORACLE_SOURCE", ROOT)) / "MarsQ_2FK"
        scratch = ROOT / "build"
        scratch.mkdir(exist_ok=True)
        cls.temporary = tempfile.TemporaryDirectory(dir=scratch, prefix="explicit-profiles-")
        cls.addClassCleanup(cls.temporary.cleanup)
        cls.work = Path(cls.temporary.name)
        source = (cls.source / "kinetic.f").read_text()
        names = ["KEQPROF", "KALPHAA", "KDIAMAG", "DFFFDPSI", "KVALIDATESPECIESNU"]
        # Include the actual reader/validation helpers, without replacing any
        # numerical routine in the executed Maxwellian path.
        names += [name for name in re.findall(r"^      SUBROUTINE +(\w+)", source, re.M)
                  if "EXPLICITPROFILE" in name or name.startswith("KPROFSPEC")
                  or name == "KREADPROFILERECORD"]
        kernel = "".join(routine(source, name) for name in dict.fromkeys(names))
        kernel += routine((cls.source / "marsq.f").read_text(), "RADF2RADP")
        metadata = "" if "KPROFSPEC" in (cls.source / "globalm.f").read_text() else "integer :: KPROFSPEC\n"
        validate = "" if "KVALIDATEEXPLICITPROFILES" not in names else """
 if(CASE_ID==3) then
  IFOWP=1;IFOWT=1
  call KVALIDATEEXPLICITPROFILES(.false.)
  IFOWP=0;IFOWT=0
 endif
 call KVALIDATEEXPLICITPROFILES(.true.)
"""
        observe = "" if "KWRITEEXPLICITPROFILES" not in names else OBSERVE
        (cls.work / "kernel.f").write_text(kernel)
        (cls.work / "driver.f90").write_text(DRIVER.replace("!METADATA", metadata)
                                             .replace("!VALIDATE", validate)
                                             .replace("!OBSERVE", observe))
        (cls.work / "mpi_fixture.f90").write_text(
            "module MPIENV\ninteger :: ISMPIRUN,RANK\ninteger,parameter :: ROOT=0\nend module\n")
        modules = ["globalm.f", "dimensim.f", "rcomdm.f", "kineticm.f",
                   "anisotropicm.f", "torquem.f", "reorbitm.f"]
        command = [os.environ.get("FC", "gfortran"), "-O0", "-fcheck=all",
                   "-fdefault-real-8", "-fdefault-double-8", "-ffixed-line-length-none",
                   "-ffree-line-length-none", "-I" + str(cls.source)]
        command += [str(cls.source / name) for name in modules]
        command += ["mpi_fixture.f90", str(cls.source / "splines.f"),
                    "kernel.f", "driver.f90", "-o", "oracle"]
        result = subprocess.run(command, cwd=cls.work, capture_output=True, text=True, timeout=60)
        if result.returncode:
            raise RuntimeError(result.stdout + result.stderr)

    def execute(self, records=None, overrides="", mode=1):
        path = self.work / "PROFSPECIES.IN"
        path.unlink(missing_ok=True)
        (self.work / "KINETIC_PROFILES.OUT").unlink(missing_ok=True)
        if records is not None:
            path.write_text("\n".join(records) + "\n")
        result = subprocess.run([str(self.work / "oracle")], cwd=self.work,
                                input=f"&ORACLE_INPUT KPROFSPEC={mode}, {overrides} /\n",
                                capture_output=True, text=True, timeout=5)
        return result

    def pipeline_rows(self, result):
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("PIPELINE_PASS", result.stdout, result.stderr)
        self.assertNotIn("STOP", result.stderr)
        rows = [list(map(float, line.split()[1:])) for line in result.stdout.splitlines()
                if line.startswith("ROW ")]
        self.assertEqual(len(rows), 3 * (2 * N - 1))
        return rows

    def close(self, actual, expected, relative=3e-12, absolute=1e-28):
        self.assertTrue(math.isclose(actual, expected, rel_tol=relative, abs_tol=absolute),
                        f"{actual:.17g} != {expected:.17g}")

    def readback_rows(self):
        path = self.work / "KINETIC_PROFILES.OUT"
        self.assertTrue(path.is_file(), "Maintained-background readback missing")
        lines = path.read_text().splitlines()
        self.assertIn("# MARS_KINETIC_PROFILES_OUT_V1", lines)
        self.assertIn("# APPROXIMATION MAINTAINED_FROZEN_FIELD", lines)
        self.assertIn("# NO_KINETIC_FORCE_BALANCE_CLAIM", lines)
        rows = [list(map(float, line.split())) for line in lines
                if line.strip() and not line.startswith("#")]
        self.assertEqual(len(rows), 3 * (2 * N - 1))
        self.assertTrue(all(len(row) == 18 for row in rows))
        return rows

    def test_flat_profiles_preserve_carrier_and_feed_native_maxwellians(self):
        rows = self.pipeline_rows(self.execute(profile()))
        for row in rows:
            species, grid, js = map(int, row[:3])
            cs, density, energy, pressure, prep, omegan, omegat, eps, cpsi, denf, prefrac, nu, omegae = row[3:]
            k = species - 1
            self.close(density, DENSITIES[k] / NREF)
            self.close(energy, NREF * ENERGIES[k] / PREF)
            self.close(pressure, DENSITIES[k] * ENERGIES[k] / PREF)
            self.close(eps, energy)
            self.close(cpsi, pressure)
            self.close(denf, DENSITIES[k] / DENSITIES[1])
            self.close(prefrac, DENSITIES[k] * ENERGIES[k] /
                       sum(n * t for n, t in zip(DENSITIES, ENERGIES)))
            for derivative in (prep, omegan, omegat):
                self.assertLess(abs(derivative), 1e-11)
            self.close(nu, (100, 200, 300)[k] * 3e-6)
            self.close(omegae, 4500 * 3e-6)
        for row in self.readback_rows():
            k = int(row[0]) - 1
            self.close(row[5], DENSITIES[k])
            self.close(row[6], ENERGIES[k])
            self.close(row[7], DENSITIES[k] * ENERGIES[k])
            for derivative in (*row[8:11], *row[13:15]):
                self.assertLess(abs(derivative), 1e-11)
            self.assertLess(abs(row[11] / DENSITIES[k]), 1e-11)
            self.assertLess(abs(row[12] / ENERGIES[k]), 1e-11)
            self.close(row[15], 4500)
            self.close(row[16], (100, 200, 300)[k])
        background = [list(map(float, line.split()[2:])) for line in
                      (self.work / "KINETIC_PROFILES.OUT").read_text().splitlines()
                      if line.startswith("# BACKGROUND ") and not line.startswith("# BACKGROUND GRID")]
        self.assertEqual(len(background), 2 * N - 1)
        total_pressure = sum(n * t for n, t in zip(DENSITIES, ENERGIES))
        for row in background:
            self.assertEqual(len(row), 11)
            cs, stor, carrier_p, kinetic_p, delta_p, carrier_dp, kinetic_dp, delta_dp, rho = row[2:]
            self.close(stor, cs**2)
            self.close(carrier_p, (.4 + .3 * stor) * PREF)
            self.close(kinetic_p, total_pressure)
            self.close(delta_p, total_pressure - carrier_p)
            self.assertGreater(abs(delta_p), 1e6)
            self.close(carrier_dp, .3)
            self.assertLess(abs(kinetic_dp), 1e-11)
            self.close(delta_dp, -.3, relative=3e-10)
            self.close(rho, 1.1 + .2 * stor)

    def test_affine_flux_density_has_analytic_native_derivatives(self):
        slope = .4
        rows = self.pipeline_rows(self.execute(profile(slope)))
        for row in rows:
            species, grid, js = map(int, row[:3])
            cs, density, energy, pressure, prep, omegan, omegat, eps, cpsi = row[3:12]
            k = species - 1
            self.close(density, DENSITIES[k] * (1 + slope * cs**2) / NREF)
            self.close(energy, NREF * ENERGIES[k] / PREF)
            self.close(prep, DENSITIES[k] * ENERGIES[k] * slope / PREF,
                       relative=3e-10)
            self.assertLess(abs(omegat), 1e-11)
            # Native logarithmic derivatives extrapolate at the integer
            # endpoints; every physical half-grid and interior whole point
            # has the exact affine-flux logarithmic derivative.
            if grid == 2 or 1 < js < N:
                expected = -slope / (1 + slope * cs**2) * 1.3 / 50 * energy / (1, -1, 1)[k]
                self.close(omegan, expected, relative=3e-10)
        for row in self.readback_rows():
            species, grid, js = map(int, row[:3])
            k = species - 1
            self.close(row[4], row[3]**2)
            self.close(row[11], DENSITIES[k] * slope, relative=3e-10)
            self.assertLess(abs(row[12]), 1e-26)
            if grid == 2 or 1 < js < N:
                self.close(row[9], slope / (1 + slope * row[4]), relative=3e-10)

    def test_numerical_drive_readback_uses_native_consumer_arrays(self):
        for label, records in (("flat", profile()), ("joint gradients", profile(.4, .6))):
            with self.subTest(background=label):
                result = self.execute(records)
                rows = self.pipeline_rows(result)
                native = {tuple(map(int, row[:3])): row for row in rows}
                consumed = {tuple(map(int, fields[:3])): tuple(map(float, fields[3:]))
                            for line in result.stdout.splitlines() if line.startswith("CONSUMED ")
                            for fields in [line.split()[1:]]}
                self.assertEqual(len(consumed), 3 * (2 * N - 1))
                readback = self.readback_rows()
                for row in readback:
                    actual = native[tuple(map(int, row[:3]))]
                    self.close(row[5], actual[4] * NREF)
                    self.close(row[6], actual[5] * PREF / NREF)
                lines = (self.work / "KINETIC_PROFILES.OUT").read_text().splitlines()
                self.assertIn("# KFASTRUN 0", [" ".join(line.split()) for line in lines])
                numerical = [list(map(float, line.split()[2:])) for line in lines
                             if line.startswith("# NUMERICAL_DRIVE ")
                             and not line.startswith("# NUMERICAL_DRIVE KP")]
                self.assertEqual(len(numerical), 3 * (2 * N - 1))
                differences = []
                raw = {tuple(map(int, row[:3])): row for row in readback}
                for row in numerical:
                    self.assertEqual(len(row), 5)
                    key = tuple(map(int, row[:3]))
                    dcd, deps = consumed[key]
                    # With constant q and DPSIDS=2CS the native-to-stor
                    # Jacobian is 1. These are the actual arrays consumed
                    # by numerical Maxwellian KI, independently printed
                    # directly from the production arrays.
                    self.close(row[3], dcd - deps, absolute=1e-14)
                    self.close(row[4], deps, absolute=1e-14)
                    differences.append(abs(row[3] - raw[key][9]))
                if label == "flat":
                    self.assertEqual(max(differences), 0)
                else:
                    # Applying a discrete derivative to n*T does not obey
                    # the analytic product rule for these nonflat inputs.
                    self.assertGreater(max(differences), 1e-6)

    def test_early_validation_ignores_stale_derived_orbit_flags(self):
        rows = self.pipeline_rows(self.execute(profile(), "CASE_ID=3"))
        for row in rows:
            self.close(row[4], DENSITIES[int(row[0]) - 1] / NREF)

    def test_invalid_files_and_unsupported_modes_reject(self):
        cases = [("missing", None)]
        mutations = {
            "version": (0, "MARS_KINETIC_PROFILES_V2"),
            "model": (1, "SELF_CONSISTENT NATIVE_CS DENSITY_M3 TEMPERATURE_ENERGY_J"),
            "units": (1, "MAINTAINED_FROZEN_FIELD NATIVE_CS DENSITY_M3 TEMPERATURE_EV"),
            "counts": (2, "2 9"),
            "duplicate identity": (4, f"1 -1 {MASSES[1]:.17e}"),
            "charge": (3, f"1 2 {MASSES[0]:.17e}"),
            "mass": (3, "1 1 1e-27"),
            "identity suffix": (3, f"1 1 {MASSES[0]:.17e} 7"),
        }
        for label, (index, replacement) in mutations.items():
            records = profile()
            records[index] = replacement
            cases.append((label, records))
        for label, token, index in (("NaN", "NaN", 1), ("Inf", "Infinity", 2),
                                    ("negative n", "-1", 1), ("zero T", "0", 2),
                                    ("nonneutral", "2e19", 1), ("wrong grid", ".01", 0),
                                    ("overflow pressure", "1e308", 2)):
            records = profile()
            fields = records[8].split()
            fields[index] = token
            records[8] = " ".join(fields)
            cases.append((label, records))
        records = profile()
        cases += [("missing row", records[:-1]), ("duplicate row", records[:9] + [records[8]] + records[9:]),
                  ("trailing record", records + ["7"]), ("trailing blank", records + [""]),
                  ("trailing field", records[:-1] + [records[-1] + " 7"]),
                  ("oversized record", records[:-1] + [records[-1] + " " * 4096])]
        for label, records in cases:
            with self.subTest(file=label):
                result = self.execute(records)
                self.assertNotIn("PIPELINE_PASS", result.stdout)
                self.assertRegex(result.stdout + result.stderr, r"(?i)PROFILE|PROFSPECIES")
        guards = ["KPERTREAD=0", "IPERTURB=0", "INCKIN=0", "KNTV=11", "KNUSPECIES=0",
                  "KDWKREAD=1", "IFOW=1", "IFOWP=1", "IFOWT=1", "IDIAMV=1", "IDIAMB=1",
                  "IDIAMTI=1", "IDIAMTE=1", "INCDPHI=1", "IGAMMADIVV=1", "NDENEQ=1",
                  "KPROFTAUTH=1", "NPROFIE=4", "NEXPV=1", "NPROFUI=1", "NPROFUE=1",
                  "SLEFT=.1", "NPROFWE=1", "ISMPIRUN=1", "CASE_ID=2", "NSPECIES=1", "NSPECIES=11", "ZTAUA0=0",
                  "NR=2,NRP1=3", "PSPECIES_FOWP(3)=1", "PSPECIES_FOWT(3)=-1"]
        for override in guards:
            with self.subTest(guard=override):
                result = self.execute(profile(), override)
                self.assertNotIn("PIPELINE_PASS", result.stdout)
                self.assertRegex(result.stdout + result.stderr, r"(?i)KPROFSPEC|EXPLICIT.*PROFILE")
        result = self.execute(profile(), mode=2)
        self.assertNotIn("PIPELINE_PASS", result.stdout)
        self.assertIn("KPROFSPEC", result.stdout + result.stderr)

    def test_default_off_uses_native_legacy_pressure_and_density_formulas(self):
        rows = self.pipeline_rows(self.execute(None, mode=0))
        self.assertFalse((self.work / "KINETIC_PROFILES.OUT").exists())
        for row in rows:
            species, grid, js = map(int, row[:3])
            cs, density, energy, pressure, prep, omegan, omegat, eps, cpsi = row[3:12]
            k = species - 1
            psi = cs**2
            density_psi = ((N - 2) / (N - 1) if grid == 1 else
                           (N - 2.5) / (N - 1)) if js == N else psi
            expected_n = (.75, 1, .25)[k] * (1.1 + .2 * density_psi)
            expected_p = (.32, .48, .2)[k] * (.4 + .3 * psi)
            self.close(density, expected_n)
            self.close(pressure, expected_p)
            self.close(cpsi, expected_p)
            self.close(prep, (.32, .48, .2)[k] * .3, relative=3e-10)
            if js < N:
                self.close(energy, expected_p / expected_n)


OBSERVE = """
 if (KPROFSPEC==1) then
  FFF=7d0; DFFF=-9d0
  call KWRITEEXPLICITPROFILES
  if (any(FFF/=7d0).or.any(DFFF/=-9d0)) error stop 'observer modified derivative workspaces'
 endif
"""

DRIVER = """
program explicit_profiles_oracle
 use rcomdm
 use dimensim
 use globalm
 use kineticm
 use mpienv
 use anisotropicm
 implicit none
 integer :: j,k,g,CASE_ID
 real*8 :: carrier(9,6),scales(6)
!METADATA
 namelist /ORACLE_INPUT/ KPROFSPEC,KPERTREAD,IPERTURB,INCKIN,KNTV,KNUSPECIES,KDWKREAD,IFOW,IFOWP,IFOWT, &
  IDIAMV,IDIAMB,IDIAMTI,IDIAMTE,INCDPHI,IGAMMADIVV,NDENEQ,KPROFTAUTH,NPROFIE,NEXPV,NPROFUI,NPROFUE, &
  SLEFT,NPROFWE,ISMPIRUN,CASE_ID,NSPECIES,ZTAUA0,NR,NRP1,PSPECIES_FOWP,PSPECIES_FOWT
 NR=8;NRP1=9;NSPECIES=3;ISMPIRUN=0;RANK=ROOT;CASE_ID=0
 KPROFSPEC=1;KPERTREAD=1;IPERTURB=1;INCKIN=1;KNTV=21;KNUSPECIES=1;KDWKREAD=0
 IFOW=0;IFOWP=0;IFOWT=0;IDIAMV=0;IDIAMB=0;IDIAMTI=0;IDIAMTE=0;INCDPHI=0;IGAMMADIVV=0
 NDENEQ=0;KPROFTAUTH=0;NPROFIE=0;NEXPV=0;NPROFUI=0;NPROFUE=0;SLEFT=0d0;NPROFWE=0
 ZTAUA0=3d-6;ZNE0=5d19;B0EXP=2.3d0;R0EXP=6.2d0;OMEGACI0=50d0;B0K=1.3d0
 PSPECIES_FOWP=0d0;PSPECIES_FOWT=0d0
 read(*,nml=ORACLE_INPUT)
 allocate(CS(9),CSM(9),Q(9),QM(9),DPSIDS(9),DPSIDSM(9),PEQ(9),PEQM(9),PPEQ(9),PPEQM(9),RHO(9),RHOM(9))
 allocate(GNUI(9),GNUIM(9),GNUE(9),GNUEM(9),FFF(9,2),DFFF(9,2),OMEGAE0(9,2))
 allocate(ESPECIES_DENF(9,2,3),ESPECIES_PREF(9,2,3),ESPECIES_DEN(9,2,3), &
          ESPECIES_PRE(9,2,3),ESPECIES_TEM(9,2,3),ESPECIES_PREP(9,2,3))
 allocate(EPSALPHA(9,2,3),EPSLONCA(9,2,3),AAK(9,2,3),CPSI(9,2,3),ZC1(9,2,3), &
          ALPHAA1(9,2,3),ALPHAA2(9,2,3),ALPHAA3(9,2,3),DCDPSIL(9,2,3),DCDPSIL2(9,2,3), &
          DEPSALPHADPSI(9,2,3),DEPSCDPSI2(9,2,3),DEPSALPHADPSI2(9,2,3), &
          OMEGASN(9,2,3),OMEGAST(9,2,3),OMEGASNA(9,2,3),OMEGASAA(9,2,3),OMEGASCA(9,2,3),NUEFF(9,2,3))
 allocate(B1U(9,1),X1U(9,1))
 allocate(ZDHMAXDPSI(9,2),ZDHMAXDPSI2(9,2),HMIN3DPSI(9,2),HMIN3DPSI2(9,2), &
          ZDHMINMAXDPSI(9,2),ZDHMINMAXDPSI2(9,2))
 ZDHMAXDPSI=0d0;ZDHMAXDPSI2=0d0;HMIN3DPSI=0d0;HMIN3DPSI2=0d0;ZDHMINMAXDPSI=0d0;ZDHMINMAXDPSI2=0d0
 do j=1,9
  CS(j)=sqrt(dble(j-1)/8d0);CSM(j)=sqrt((dble(j)-.5d0)/8d0)
 enddo
 DPSIDS=2d0*CS;DPSIDSM=2d0*CSM;Q=1.7d0;QM=1.7d0
 PEQ=.4d0+.3d0*CS**2;PEQM=.4d0+.3d0*CSM**2;PPEQ=.3d0;PPEQM=.3d0
 RHO=1.1d0+.2d0*CS**2;RHOM=1.1d0+.2d0*CSM**2
 carrier(:,1)=PEQ;carrier(:,2)=PEQM;carrier(:,3)=RHO;carrier(:,4)=RHOM;carrier(:,5)=PPEQ;carrier(:,6)=PPEQM
 scales=[ZNE0,B0EXP,R0EXP,OMEGACI0,ZTAUA0,B0K]
 B1U=cmplx(3d0,-4d0,8);X1U=cmplx(-2d0,5d0,8)
 ESPECIES_M=2d0;ESPECIES_M(2)=1d0/1836.15267343d0;ESPECIES_M(3)=3d0
 ESPECIES_Z=1d0;ESPECIES_Z(2)=-1d0;ISPECIES_F0=0
 if(CASE_ID==2) ISPECIES_F0(3)=1
 ESPECIES_DENF=.25d0;ESPECIES_PREF=.25d0;ESPECIES_DEN=0d0;ESPECIES_PRE=0d0;ESPECIES_TEM=0d0
 FFF=0d0;DFFF=0d0;ALPHAP=.4d0;NKSMOOTHR=0;NUMODEL=0;INUTYPE=0;ASPCT=3d0;KFASTRUN=0
 GNUI=.01d0;GNUIM=.015d0;GNUE=.1d0;GNUEM=.15d0;NUEFFSA(1:3)=[100d0,200d0,300d0]
 OMEGAE0=4500d0*ZTAUA0
!VALIDATE
 call KEQPROF
 call KALPHAA
 call KDIAMAG
!OBSERVE
 if(any(PEQ/=carrier(:,1)).or.any(PEQM/=carrier(:,2)).or.any(RHO/=carrier(:,3)).or.any(RHOM/=carrier(:,4))) &
  error stop 'carrier profiles modified'
 if(any(PPEQ/=carrier(:,5)).or.any(PPEQM/=carrier(:,6))) error stop 'carrier derivatives modified'
 if(any([ZNE0,B0EXP,R0EXP,OMEGACI0,ZTAUA0,B0K]/=scales)) error stop 'reference scales modified'
 if(any(B1U/=cmplx(3d0,-4d0,8)).or.any(X1U/=cmplx(-2d0,5d0,8))) error stop 'frozen field modified'
 do k=1,3
  do g=1,2
   do j=1,9
    if(g==2.and.j==9) cycle
    write(*,'(A,3I4,13ES25.16)') 'ROW ',k,g,j,merge(CS(j),CSM(j),g==1),ESPECIES_DEN(j,g,k), &
     ESPECIES_TEM(j,g,k),ESPECIES_PRE(j,g,k),ESPECIES_PREP(j,g,k),OMEGASN(j,g,k),OMEGAST(j,g,k), &
     EPSALPHA(j,g,k),CPSI(j,g,k),ESPECIES_DENF(j,g,k),ESPECIES_PREF(j,g,k),NUEFF(j,g,k),OMEGAE0(j,g)
    write(*,'(A,3I4,2ES25.16)') 'CONSUMED ',k,g,j,DCDPSIL(j,g,k),DEPSALPHADPSI(j,g,k)
   enddo
  enddo
 enddo
 print '(A)', 'PIPELINE_PASS'
end program
subroutine KJP_ISOEXTRA(k)
 integer,intent(in)::k
 error stop 'unexpected non-Maxwellian helper'
end subroutine
subroutine KALPHAA_TYPE3(k)
 integer,intent(in)::k
 error stop 'unexpected non-Maxwellian helper'
end subroutine
subroutine KCPSI_TYPE3(k)
 integer,intent(in)::k
 error stop 'unexpected non-Maxwellian helper'
end subroutine
subroutine KCPSI_TYPE4(k)
 integer,intent(in)::k
 error stop 'unexpected non-Maxwellian helper'
end subroutine
subroutine KCPSI_TYPE5(k)
 integer,intent(in)::k
 error stop 'unexpected non-Maxwellian helper'
end subroutine
subroutine KCPSI_TYPE6(k)
 integer,intent(in)::k
 error stop 'unexpected non-Maxwellian helper'
end subroutine
"""


if __name__ == "__main__":
    unittest.main()
