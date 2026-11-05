"""Figure of PR #131: EXAMPLE axis-grid vacuum round trip of master 0668492, the fix and master before #103.

usage: MPLBACKEND=Agg python make_131.py DATA 131.png, with DATA/<label>/{BPLASMA.OUT,result.json}
from make_131_data.py run in tests/ with the fix (label fix), 0668492 (m0668492) and de75a66 (pre103) builds.
"""
import json, sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

d = Path(sys.argv[1]); out = Path(sys.argv[2])
OI = dict(black='#000000', orange='#E69F00', sky='#56B4E9', green='#009E73', blue='#0072B2',
          verm='#D55E00', purple='#CC79A7')
VAR = [('m0668492', 'master 0668492 (before)', OI['verm'], 'o', '-'),
       ('fix', 'this PR', OI['blue'], 's', '-'),
       ('pre103', 'master before #103 (de75a66)', OI['orange'], 'x', '--')]


def load(label):
    p = d / label / 'BPLASMA.OUT'
    rows = p.read_text().split('\n')
    nm, ns = (int(float(x)) for x in rows[0].split()[:2])
    modes = [int(round(float(rows[1 + i].split()[0]))) for i in range(nm)]
    body = np.loadtxt(p, skiprows=1 + nm)
    b = (body[:, 0] + 1j * body[:, 1]).reshape(nm, ns).T
    return dict(zip(modes, b.T)), json.loads((d / label / 'result.json').read_text())


data = {v[0]: load(v[0]) for v in VAR}
res0 = data['fix'][1]
ibnm = res0['npsi'] + res0['kv_bnm'] - 1   # 0-based row of IBNM
npl = res0['npsi']                          # plasma boundary row
presc = {int(m): complex(*v) for m, v in res0['prescribed'].items()}
fig, (a, b) = plt.subplots(1, 2, figsize=(10.5, 4.2))
ms = sorted(int(m) for m in res0['b1_ibnm'])
for lab, name, col, mk, ls in VAR:
    r = data[lab][1]
    err = [abs(complex(*r['b1_ibnm'][str(m)]) - presc.get(m, 0j)) for m in ms]
    a.semilogy(ms, err, color=col, marker=mk, ls=ls, mfc='none' if mk != 'x' else None,
               label=f"error, {name}: misfit {r['misfit']:.4f}")
a.semilogy(sorted(presc), [abs(presc[m]) for m in sorted(presc)], 'k*', ms=11, ls='none',
           label=r'prescribed $|b_{1,m}|$ (BS_ESC input)')
a.set_xlabel('poloidal harmonic $m$')
a.set_ylabel(r'$|b_{1,m}|$ or $|b_{1,m}^{\rm vac} - b_{1,m}^{\rm presc}|$ on IBNM (MARS normalised units)')
a.set_title('(a) normal field returned on IBNM, $n=1$'); a.set_ylim(2e-8, 1e-1)
a.legend(fontsize=7.5, loc='upper left'); a.grid(alpha=.3, which='both')
fix = data['fix'][0]
cols = [OI['verm'], OI['green'], OI['purple'], OI['sky']]
mks = ['o', 's', '^', 'D']
for (m, col, mk) in zip((0, 1, 2, 3), cols, mks):
    y = np.abs(data['m0668492'][0][m] - fix[m])[1:ibnm + 6]
    k = np.arange(1, ibnm + 6)
    b.semilogy(k, y, color=col, marker=mk, markevery=4, mfc='none', label=f'$m={m}$')
same = max(np.abs(data['pre103'][0][m] - fix[m]).max() for m in fix)
b.axvline(npl, color='0.5', lw=.8, ls=':')
b.axvline(ibnm, color='0.5', lw=.8, ls=':')
b.text(npl - .5, 8e-6, 'plasma boundary', fontsize=7, ha='right')
b.text(ibnm + .5, 8e-6, 'IBNM', fontsize=7)
b.set_ylim(1e-10, 1e-4)
b.set_xlabel('radial grid surface index (0 = magnetic axis)')
b.set_ylabel(r'$|b_{1,m}^{\rm 0668492} - b_{1,m}^{\rm PR}|$ (MARS normalised units)')
b.set_title('(b) vacuum field before minus with this PR')
b.text(0.02, 0.03, ('master before #103 and this PR: identical, all $m$ and surfaces' if same == 0 else f'master before #103 vs this PR: max difference {same:.1e}'),
       transform=b.transAxes, fontsize=7.5)
b.legend(fontsize=8, loc='upper left', ncol=4); b.grid(alpha=.3, which='both')
fig.suptitle('INCFEED=4 vacuum solve of the ESC currents fitted to a prescribed normal field (EXAMPLE equilibrium)', fontsize=10)
fig.tight_layout(); fig.savefig(out, dpi=150)
for lab, *_ in VAR:
    print(lab, data[lab][1]['misfit'])
