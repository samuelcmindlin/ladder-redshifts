"""Batch mode for real ALMA spectra: every source with cube cut-outs under data/ is turned into an
aperture spectrum, run through the engine and scored against data/truth.csv.

Discovery: any data/**/<TARGET>_spwNN.npz or <TARGET>_fNN.NNN.npz (written by
scripts/alma/extract_cube_spectra.py). Cut-outs centred on a known source position (extractor --targets) give a
point-source spectrum at the centre; others get the detection-grown aperture of realspec.py. Either way the
per-channel noise is put on the scale of the spectrum's own line-free channels (realspec.calibrate_noise; the
factors go to data/noise_calibration.csv). If a source's cut-outs are absent (they are not shipped in the
repo), its saved spectrum data/spectra/<TARGET>.csv is used instead. Spectra built from ADMIT summaries
(scripts/alma/admit_spectra.py) carry an "admit" column: windows measured at the source vs zero-filled
non-detections. For those, the engine is also run on the measured windows alone, to show what the
non-detections contribute.

truth.csv: source, z_published, prior, reference, notes, status ("confirmed" or "tentative"; a tentative
value is a published but unconfirmed line identification, and verdicts against it are marked as such).

Writes data/spectra/<TARGET>.csv, results/redshift_real_batch.json, results/redshift_real_batch.csv and
figures/redshift_real_batch.png (one row per source: spectrum with the engine's lines + posterior candidates).
--admit scores the ADMIT-route spectra in data/admit/spectra/ instead (no cut-outs), writing the same files
with the suffix _admit.
Usage: python -m ladder_redshifts.batch_real [--prior auto|uniform|spt|herschel]
"""
from .paths import AUDIO, FIG, RESULTS, DATA  # noqa: E402
import argparse
import csv
import json
import re
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .ladder import Grid, line_list, prior_spt, prior_herschel, TOL_LN, CONF_WIN  # noqa: E402
from .realspec import (aperture, spectrum, spectrum_point, run_engine, grid_windows, noise_scale,  # noqa: E402
                       record_noise_scales)

INK, INK2, GRID_C, SURF = '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'
COL = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300']
PRIORS = {'uniform': None, 'spt': prior_spt, 'herschel': prior_herschel}


def discover():
    groups = defaultdict(list)
    for f in sorted((DATA).rglob('*.npz')):
        m = re.match(r'(.+)_(spw\d+|f\d+\.\d+)$', f.stem)
        if m and 'collapsed' not in f.stem:
            groups[m.group(1)].append(f)
    for f in sorted((DATA / 'spectra').glob('*.csv')) if (DATA / 'spectra').exists() else []:
        groups.setdefault(f.stem, [])
    return groups


def load_files(files):
    spws = []
    for f in files:
        a = np.load(f, allow_pickle=False)
        hdr = json.loads(str(a['header'])) if 'header' in a.files else {}
        o = np.argsort(a['freq_hz'])                      # lower-sideband windows come in descending order
        spws.append(dict(name=f.stem, freq=a['freq_hz'][o] / 1e9, cut=a['cutout'][o].astype(np.float64),
                         rms=a['rms_jyb'][o].astype(np.float64), beams=a['beams'], pix=float(a['pix_arcsec']),
                         k=int(a['chan_bin']), centred=hdr.get('centre_from') == 'target position'))
    spws.sort(key=lambda s: s['freq'].min())
    return spws


def read_spectrum(path):
    with open(path, newline='') as fh:
        rows = list(csv.DictReader(fh))
    out = []
    for name in dict.fromkeys(r['spw'] for r in rows):
        rr = [r for r in rows if r['spw'] == name]
        out.append(dict(spw=name, freq=np.array([float(r['freq_ghz']) for r in rr]),
                        S=np.array([float(r['flux_mjy']) for r in rr]),
                        sig=np.array([float(r['noise_mjy']) for r in rr]),
                        cont=np.array([float(r['continuum_mjy']) for r in rr]),
                        kind=rr[0].get('admit') or 'measured'))
    return out


def write_spectrum(spec, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = ['spw,freq_ghz,flux_mjy,noise_mjy,continuum_mjy']
    for w in spec:
        for f, s, n, c in zip(w['freq'], w['S'], w['sig'], w['cont']):
            rows.append(f"{w['spw']},{f:.6f},{s:.5f},{n:.5f},{c:.5f}")
    path.write_text('\n'.join(rows) + '\n')


def to_grid_auto(spec):
    """Grid spanning the spectrum's coverage (3% margin); windows combined order-independently
    (realspec.grid_windows)."""
    lo = min(np.nanmin(w['freq']) for w in spec) * 0.97
    hi = max(np.nanmax(w['freq']) for w in spec) * 1.03
    return grid_windows(spec, Grid(lo, hi))


def read_truth():
    path = DATA / 'truth.csv'
    if not path.exists():
        return {}
    with open(path, newline='') as fh:
        return {r['source']: r for r in csv.DictReader(fh)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--prior', default='auto', choices=['auto', 'uniform', 'spt', 'herschel'])
    ap.add_argument('--admit', action='store_true', help='score data/admit/spectra/*.csv (ADMIT route)')
    a = ap.parse_args()
    truth = read_truth()
    sdir = DATA / ('admit/spectra' if a.admit else 'spectra')
    sfx = '_admit' if a.admit else ''
    groups = {f.stem: [] for f in sorted(sdir.glob('*.csv'))} if a.admit else discover()
    rows, plots = [], []
    for src, files in sorted(groups.items()):
        csvp = sdir / f'{src}.csv'
        if files:
            spws = load_files(files)
            if all(s['centred'] for s in spws):         # cut-outs centred on the source's known position
                spec = spectrum_point(spws)
                route = f'cut-outs ({len(files)} windows, point source at the published position)'
            else:
                img, csig, mask, circ = aperture(spws)
                spec = spectrum(spws, mask)
                route = f'cut-outs ({len(files)} windows, aperture {int(mask.sum())} px)'
            write_spectrum(spec, csvp)
            record_noise_scales(csvp.relative_to(DATA).as_posix(), spec)    # factors applied by calibrate_noise
        elif csvp.exists():
            spec = read_spectrum(csvp); route = 'saved spectrum'
            kinds = [w.get('kind', 'measured') for w in spec]
            if any(k != 'measured' for k in kinds):
                route = (f"ADMIT summary ({kinds.count('measured')} windows measured at the source, "
                         f"{kinds.count('zero-filled')} zero-filled non-detections)")
        else:
            continue
        scat = [noise_scale(w['S'], w['sig']) for w in spec if w.get('kind', 'measured') == 'measured']
        scat = float(np.nanmedian(scat)) if np.isfinite(scat).any() else None
        if scat is not None and abs(scat - 1) > 0.15:
            print(f"{src:14s} note: line-free flux/noise scatters by {scat:.2f}, not 1; the noise column looks "
                  f"mis-scaled (realspec.calibrate_noise corrects it)")
        t = truth.get(src, {})
        pname = a.prior if a.prior != 'auto' else (t.get('prior') or ('spt' if src.startswith('SPT') else 'uniform'))
        zt = float(t['z_published']) if t.get('z_published') else None
        if not any(w.get('kind', 'measured') == 'measured' for w in spec):
            rows.append(dict(source=src, route=route, prior=pname, z_published=zt, z_engine=None, confidence=None,
                             verdict='no spectrum on the source (ADMIT)', z_flat=None, conf_flat=None, z_sudoku=None,
                             z_naive=None, truth_status=t.get('status') if zt is not None else None, p_published=None,
                             reference=t.get('reference', '')))
            print(f"{src:14s} no spectrum on the source: ADMIT's spectra all fall elsewhere in the field "
                  f"(needs cubes or cut-outs); published {zt}")
            continue
        r = run_engine(*to_grid_auto(spec), prior=PRIORS[pname])
        meas = [w for w in spec if w.get('kind', 'measured') == 'measured']
        r_meas = run_engine(*to_grid_auto(meas), prior=PRIORS[pname]) if 0 < len(meas) < len(spec) else None
        verdict = 'no truth'
        if zt is not None:
            ok = abs(np.log1p(r['sled'][0]) - np.log1p(zt)) < TOL_LN
            verdict = ('confident & correct' if ok else 'CONFIDENTLY WRONG') if r['sled'][1] >= 0.95 else \
                      ('ambiguous (best = truth)' if ok else 'ambiguous (best != truth)')
            if t.get('status', 'confirmed') == 'tentative':
                verdict += ' [tentative ID]'
        cov = [round(min(w['freq'].min() for w in spec), 2), round(max(w['freq'].max() for w in spec), 2)]
        row = dict(source=src, route=route, coverage_ghz=cov, prior=pname, z_published=zt,
                   z_engine=round(r['sled'][0], 4), confidence=round(r['sled'][1], 3), verdict=verdict,
                   z_flat=round(r['flat'][0], 4), conf_flat=round(r['flat'][1], 3),
                   z_sudoku=None if r['sudoku'][0] is None else round(r['sudoku'][0], 4), n_sudoku=r['sudoku'][1],
                   z_naive=None if r['naive'] is None else round(r['naive'], 3),
                   detections=[(round(f, 3), round(s, 1)) for f, s in r['detections']],
                   candidates=[(round(z, 4), round(m, 3)) for z, m in r['candidates'][:5]],
                   truth_status=t.get('status', 'confirmed') if zt is not None else None,
                   p_published=None if zt is None else
                   round(float(r['post'][np.abs(r['lnz'] - np.log1p(zt)) < CONF_WIN].sum()), 3),
                   reference=t.get('reference', ''))
        if r_meas is not None:
            row['measured_windows_only'] = dict(
                z_engine=round(r_meas['sled'][0], 4), confidence=round(r_meas['sled'][1], 3),
                candidates=[(round(z, 4), round(m, 3)) for z, m in r_meas['candidates'][:6]])
        rows.append(row); plots.append((src, spec, r, zt))
        print(f"{src:14s} {verdict:26s} engine z = {row['z_engine']:.4f} (conf {row['confidence']:.3f}); "
              f"published {zt}; Sudoku {row['z_sudoku']}; naive {row['z_naive']}; [{route}; prior {pname}]")
        print(f"{'':14s} candidates {row['candidates'][:4]}")
        if r_meas is not None:
            mo = row['measured_windows_only']
            print(f"{'':14s} measured windows only: z = {mo['z_engine']:.4f} (conf {mo['confidence']:.3f}); "
                  f"candidates {mo['candidates'][:4]}")
    if not rows:
        print('no real sources found under data/'); return
    (RESULTS / f'redshift_real_batch{sfx}.json').write_text(json.dumps(rows, indent=1, default=float))
    keys = ['source', 'z_published', 'truth_status', 'z_engine', 'confidence', 'p_published', 'verdict', 'z_flat',
            'conf_flat', 'z_sudoku', 'z_naive', 'prior', 'reference']
    with open(RESULTS / f'redshift_real_batch{sfx}.csv', 'w', newline='') as fh:
        wr = csv.writer(fh)
        wr.writerow(keys)
        for row in rows:
            wr.writerow(['' if row.get(k) is None else row[k] for k in keys])
    summ = {v: sum(r['verdict'] == v for r in rows) for v in dict.fromkeys(r['verdict'] for r in rows)}
    print('\nsummary:', summ)
    amb = [r for r in rows if str(r['verdict']).startswith('ambiguous')]
    if amb:
        print(f"ambiguous answers: mean stated confidence {np.mean([r['confidence'] for r in amb]):.2f}; best candidate "
              f"= published z in {sum('best = truth' in r['verdict'] for r in amb)} of {len(amb)}; mean posterior on "
              f"the published z {np.mean([r['p_published'] for r in amb]):.2f}")
    plots = [q for q in plots if q[2] is not None]
    n = len(plots)
    fig, axs = plt.subplots(n, 2, figsize=(15, 2.6 * n + 0.6), facecolor=SURF, squeeze=False,
                            gridspec_kw=dict(width_ratios=[3, 1]))
    names, nus = line_list()
    for i, (src, spec, r, zt) in enumerate(plots):
        ax = axs[i, 0]
        for w in spec:
            if w.get('kind', 'measured') == 'measured':
                ax.step(w['freq'], w['S'] / w['sig'], where='mid', color=COL[0], lw=0.7)
            else:                                   # zero-filled ADMIT non-detection: show coverage only
                ax.axvspan(w['freq'].min(), w['freq'].max(), color='#ebeae6', lw=0, zorder=0)
        zb = r['sled'][0]
        lo = min(w['freq'].min() for w in spec); hi = max(w['freq'].max() for w in spec)
        for nm, nu in zip(names, nus):
            f = nu / (1 + zb)
            if lo - 0.3 < f < hi + 0.3:
                ax.axvline(f, color=COL[1], lw=0.9, ls='--')
                ax.text(f, 0.97, nm, rotation=90, fontsize=6.5, color=COL[1], ha='right', va='top',
                        transform=ax.get_xaxis_transform())
        ax.set_xlim(lo - 0.3, hi + 0.3)
        ax.set_ylabel('S/N', fontsize=8, color=INK2)
        tent = ' (tentative)' if truth.get(src, {}).get('status') == 'tentative' else ''
        extra = '; grey = window covered, no line found by ADMIT' if any(w.get('kind', 'measured') != 'measured' for w in spec) else ''
        ax.set_title(f'{src}: engine z = {zb:.4f} (confidence {r["sled"][1]:.3f}); published {zt}{tent}{extra}',
                     fontsize=9, color=INK, loc='left')
        ax2 = axs[i, 1]
        cands = r['candidates'][:5]
        for j, (z, m) in enumerate(cands):
            good = zt is not None and abs(np.log1p(z) - np.log1p(zt)) < CONF_WIN
            ax2.barh(j, m, color=COL[2] if good else '#c9c8c3', height=0.6)
            ax2.text(m + 0.02, j, f'z={z:.3f}', va='center', fontsize=7.5, color=INK2)
        ax2.set_xlim(0, 1.3); ax2.set_yticks([]); ax2.invert_yaxis()
        ax2.set_title('posterior mass per candidate (green = published z)', fontsize=8, color=INK, loc='left')
        for a_ in (ax, ax2):
            a_.set_facecolor(SURF)
            for sp in ('top', 'right'):
                a_.spines[sp].set_visible(False)
            a_.tick_params(colors=INK2, labelsize=7.5)
    axs[-1, 0].set_xlabel('observed frequency (GHz)', fontsize=8, color=INK2)
    fig.tight_layout()
    fig.savefig(FIG / f'redshift_real_batch{sfx}.png', dpi=120, facecolor=SURF)


if __name__ == '__main__':
    main()
