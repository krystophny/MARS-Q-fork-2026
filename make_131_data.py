"""EXAMPLE axis-grid vacuum round trip (tests/test_axis_vacuum_roundtrip.py) for one MARS executable.

usage (in tests/, PYTHONPATH=.): python make_131_data.py LABEL MARS_EXE OUTDIR  -> OUTDIR/LABEL/{BPLASMA.OUT,result.json}
"""
import json, os, shutil, sys, tempfile
from pathlib import Path

label, exe, out = sys.argv[1], sys.argv[2], Path(sys.argv[3])
os.environ['MARS_EXE'] = exe
import test_axis_vacuum_roundtrip as t

case = t.AxisVacuumRoundtrip('test_vacuum_solve_reproduces_the_fitted_normal_field')
with tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR')) as tmp:
    case.work = Path(tmp)
    misfit = case.roundtrip_misfit(t.MMAX)
    dest = out / label
    dest.mkdir(parents=True, exist_ok=True)
    for name in ('BPLASMA.OUT', 'RUN.IN'):
        shutil.copy2(case.work / 'vac' / name, dest / name)
    shutil.copy2(case.work / 'esc' / 'INPUT_BNM.IN', dest / 'INPUT_BNM.IN')
    got = t.b1_on_surface(case.work / 'vac' / 'BPLASMA.OUT', t.NPSI + t.KV_BNM - 1)
res = dict(label=label, exe=exe, misfit=misfit, npsi=t.NPSI, kv_bnm=t.KV_BNM, mmax=t.MMAX,
           prescribed={str(m): [complex(v).real, complex(v).imag] for m, v in t.BNM.items()},
           b1_ibnm={str(m): [v.real, v.imag] for m, v in got.items()})
(dest / 'result.json').write_text(json.dumps(res, indent=1))
print(label, 'misfit', misfit)
