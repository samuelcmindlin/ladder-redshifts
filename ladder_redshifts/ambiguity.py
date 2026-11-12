"""Why ladder redshifts fail by octaves: the candidate fundamentals form an undertone series.

Given two lines with frequency ratio b/a (lowest terms), the CO assignments (J_a, J_b) = (k a, k b)
all fit, for k = 1, 2, 3, ...: candidate fundamentals F/k (an undertone series), candidate
redshifts 1+z_k = k (1+z_1) (an overtone series in 1+z). What breaks the tie is (i) predicted
rungs that fall inside the coverage but are absent (the 'undertone check'), (ii) off-ladder
lines ([CI], [CII]), (iii) the SLED 'timbre', (iv) - in principle - the ladder's tiny anharmonicity.

Panels: A, the HerBS-108 octave family on a log-frequency axis; B, the ladder's ambiguity function
(expected likelihood vs trial redshift ratio) for two setups; C, the anharmonicity tie-breaker.
Also writes audio/herbs108_octave_ambiguity.wav: the same ambiguity, heard as pitch.
"""
from .paths import AUDIO, FIG, RESULTS, DATA  # noqa: E402
import json
import numpy as np
from scipy.io import wavfile
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .ladder import (CO, CI10, CI21, C_KMS, Grid, Engine, line_list, sled_flux, simulate,  # noqa: E402
                    matched_filter, DLN)

INK, INK2, GRID_C, SURF = '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'
COL = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300']

OBS = (97.75, 146.60)          # HerBS-108 (z-GAL); 146.60 = CO(6-5) at z = 3.7166
COVER = [(75.874, 106.868), (127.377, 158.380)]
Z_TRUE = 3.7166


def family(nu_a, nu_b, kmax=3):
    r = nu_b / nu_a
    a, b = next((a, b) for b in range(2, 14) for a in range(1, b) if abs(b / a - r) < 1e-3 and np.gcd(a, b) == 1)
    rows = []
    for k in range(1, kmax + 1):
        Ja, Jb = k * a, k * b
        z = np.mean([CO[Ja - 1] / nu_a, CO[Jb - 1] / nu_b]) - 1
        rows.append(dict(k=k, J=(Ja, Jb), z=float(z), f0_GHz=float(CO[0] / (1 + z))))
    return (a, b), rows


def tie_breaker():
    """Velocity mismatch between the two lines when a pair is assigned to the wrong octave."""
    out = []
    for (Ja, Jb) in ((2, 3), (1, 3), (3, 4), (2, 5), (3, 5), (4, 5)):
        for k in (2, 3):
            if k * Jb > 13:
                continue
            d = np.log(CO[k * Ja - 1] / CO[k * Jb - 1]) - np.log(CO[Ja - 1] / CO[Jb - 1])
            out.append(dict(pair=(Ja, Jb), alt=(k * Ja, k * Jb), mismatch_kms=float(abs(d) * C_KMS)))
    return out


def ambiguity(setup, zs=np.linspace(2.0, 5.0, 31), Jc=3.7):
    grid = Grid(30.0, 480.0)
    cov = grid.mask(setup)
    T = sled_flux(Jc, 0.29)
    # the ambiguity function is a property of the line set, so the lines sit exactly where predicted (dv_max = 0)
    eng = Engine(grid, cov, T[None, :], zmin=0.3, zmax=40.0, dv_max=0)
    flat = Engine(grid, cov, np.ones((1, len(T))), zmin=0.3, zmax=40.0, dv_max=0)
    rr = np.arange(np.log(0.25), np.log(4), DLN)
    acc = {'sled': [], 'flat': []}
    for z in zs:
        Y = eng.filter(noiseless(grid, z, 10 * T / T[3]))       # CO(4-3) at S/N 10, no noise; observed pixels only
        for name, e in (('sled', eng), ('flat', flat)):
            L = e.loglike(Y)
            acc[name].append(np.interp(rr, e.lnz - np.log1p(z), L / L.max(), left=np.nan, right=np.nan))
    # display: running max over the 'correct' tolerance (+/- 0.003 in ln(1+z)) so narrow peaks stay visible
    from scipy.ndimage import maximum_filter1d
    w = int(2 * 0.003 / DLN) + 1
    return np.exp(rr), {k: maximum_filter1d(np.nan_to_num(np.nanmean(np.array(v), axis=0)), w) for k, v in acc.items()}


def noiseless(grid, z, amp):
    class _Zero:
        def normal(self, mu, s, n):
            return np.zeros(n)
    return simulate(grid, z, amp, 450, _Zero())


def plot(fam, path):
    fig = plt.figure(figsize=(15, 8.2), facecolor=SURF)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1], width_ratios=[1.6, 1])
    ax = fig.add_subplot(gs[0, :])
    ax.set_facecolor(SURF)
    for lo, hi in COVER:
        ax.axvspan(lo, hi, color=GRID_C, alpha=0.8, lw=0)
    ax.axvspan(70.689, 78.433, color=COL[2], alpha=0.12, lw=0)
    ax.text(74.5, 3.75, 'extra setup\n(added later)', fontsize=7.5, color=INK2, ha='center')
    for nu in OBS:
        ax.axvline(nu, color=INK, lw=2.2)
    ax.text(OBS[0] - 0.6, 3.75, 'observed\n97.75 GHz', fontsize=8, ha='right', color=INK)
    ax.text(OBS[1] - 0.8, 3.75, 'observed\n146.60 GHz', fontsize=8, ha='right', color=INK)
    names, nus = line_list()
    for row in fam:
        y = 3 - row['k']
        zz = row['z']
        done = []
        for n, nu in zip(names, nus):
            f = nu / (1 + zz)
            if any(abs(f - g) < 1.0 for g in done):
                continue
            partner = [m for m, mu in zip(names, nus) if m != n and abs(mu / (1 + zz) - f) < 1.0]
            if partner:
                n = n + '+' + partner[0]; done.append(f)
            if 60 < f < 175:
                inside = any(lo <= f <= hi for lo, hi in COVER)
                is_obs = any(abs(f - o) < 0.3 for o in OBS)
                ci = n.startswith('[CI]')
                c = INK if is_obs else (COL[1] if inside else INK2)
                ax.plot([f, f], [y - 0.32, y + 0.32], color=c, lw=2.4 if is_obs else 1.6, ls=':' if ci else '-')
                if inside and not is_obs:
                    ax.text(f, y + 0.38, n + '\nshould be seen', fontsize=6.8, ha='center', color=COL[1])
                elif not is_obs:
                    ax.text(f, y - 0.5, n, fontsize=6.5, ha='center', color=INK2)
                else:
                    ax.text(f + 0.6, y - 0.05, n, fontsize=7, color=INK)
    ax.set_yticks([2, 1, 0])
    ax.set_yticklabels([f"k = {r['k']}: CO J = {r['J']}\nz = {r['z']:.3f}, fundamental\n{r['f0_GHz']:.2f} GHz = F/{r['k']}"
                        for r in fam], fontsize=8, color=INK)
    ax.set_xscale('log'); ax.set_xlim(62, 175); ax.set_ylim(-0.9, 4.3)
    from matplotlib.ticker import NullFormatter
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xticks([70, 80, 90, 100, 120, 140, 160]); ax.set_xticklabels(['70', '80', '90', '100', '120', '140', '160'])
    ax.tick_params(axis='y', length=0)
    ax.set_xlabel('observed frequency (GHz, log axis); grey = NOEMA coverage; orange = predicted inside coverage but not seen',
                  color=INK2, fontsize=9)
    for s in ('top', 'right', 'left'):
        ax.spines[s].set_visible(False)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.set_title('A. One pair of lines, an undertone family of answers (HerBS-108): candidate fundamentals F, F/2, F/3; '
                 'redshifts 1+z = 2.358 x (1, 2, 3).  True: k = 2.', fontsize=10, color=INK, loc='left')

    ax = fig.add_subplot(gs[1, 0]); ax.set_facecolor(SURF)
    for j, (lab, setup) in enumerate((('Band 3 scan', [(84.2, 114.9)]), ('Band 2 full (WSU)', [(67.0, 116.0)]),
                                      ('Band 3 + Band 4', [(84.2, 114.9), (125.0, 163.0)]))):
        r, acc = ambiguity(setup)
        ax.plot(r, acc['sled'], color=COL[j], lw=1.8, label=lab)
    for rr_, lab in ((2, '2'), (1.5, '3/2'), (4 / 3, '4/3'), (0.5, '1/2'), (2 / 3, '2/3'), (0.75, '3/4'), (3, '3'), (1 / 3, '1/3')):
        ax.axvline(rr_, color=INK2, lw=0.6, ls=':')
        ax.text(rr_, 1.04, lab, fontsize=7, color=INK2, ha='center')
    ax.set_xscale('log'); ax.set_xlim(0.25, 4); ax.set_ylim(-0.02, 1.1)
    ax.set_xticks([0.25, 0.5, 1, 2, 4]); ax.set_xticklabels(['1/4', '1/2', '1', '2', '4'])
    from matplotlib.ticker import NullFormatter
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel('trial (1+z) / true (1+z)', color=INK2, fontsize=9)
    ax.set_ylabel('expected log-likelihood (true = 1)', color=INK2, fontsize=9)
    ax.set_title('B. Ambiguity function of the CO ladder (noiseless, averaged over z = 2-5)', fontsize=10, color=INK, loc='left')
    ax.legend(fontsize=8, frameon=False)

    ax = fig.add_subplot(gs[1, 1]); ax.set_facecolor(SURF)
    tb = tie_breaker()
    labs = [f"{t['pair']}→{t['alt']}" for t in tb]
    ax.barh(range(len(tb)), [t['mismatch_kms'] for t in tb], color=COL[0], height=0.6)
    for i, t in enumerate(tb):
        ax.text(t['mismatch_kms'] + 1, i, f"{t['mismatch_kms']:.0f}", va='center', fontsize=7.5, color=INK2)
    ax.set_yticks(range(len(tb))); ax.set_yticklabels(labs, fontsize=7.5)
    snr = np.array([10, 30, 100])
    for s_, ls in zip(snr, (':', '--', '-')):
        sig = np.sqrt(2) * (np.sqrt(2) * 450 / 2.3548 / s_)
        ax.axvline(sig, color=COL[1], lw=1.4, ls=ls, label=f'1σ on the pair difference, both lines at S/N {s_}')
    ax.invert_yaxis(); ax.set_xlim(0, 190)
    ax.legend(fontsize=7, frameon=False, loc='lower right')
    ax.set_xlabel('velocity mismatch between the two lines (km/s)', color=INK2, fontsize=9)
    ax.set_title('C. Anharmonic tie-breaker: wrong harmonic family vs centroid precision', fontsize=10, color=INK, loc='left')
    for a in fig.axes[1:]:
        for s in ('top', 'right'):
            a.spines[s].set_visible(False)
        a.tick_params(colors=INK2, labelsize=8); a.grid(color=GRID_C, lw=0.6); a.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(path, dpi=130, facecolor=SURF)


def sonify(path, sr=44100):
    """Map GHz -> Hz so the HerBS-108 lines become 440 and 660 Hz. Two lines alone: the ear hears a
    virtual pitch of 220 Hz (k = 1). Add the decisive third line (73.3 GHz -> 330 Hz) and the
    virtual pitch drops an octave to 110 Hz (k = 2) - exactly how the redshift was resolved."""
    scale = 440.0 / OBS[0]
    def tone(freqs, amps, dur, fade=0.05):
        t = np.arange(int(dur * sr)) / sr
        y = sum(a * np.sin(2 * np.pi * f * t) for f, a in zip(freqs, amps))
        e = np.minimum(1, np.minimum(t / fade, (dur - t) / fade))
        return y * e
    gap = np.zeros(int(0.5 * sr))
    obs = [OBS[0] * scale, OBS[1] * scale]                       # 440, 660 Hz
    third = 73.315 * scale                                       # 330 Hz
    f1 = CO[0] / (1 + 1.3583) * scale                            # 220 Hz
    f2 = CO[0] / (1 + Z_TRUE) * scale                            # 110 Hz
    seq = [tone(obs, [1, 1], 3.0), gap,
           tone([f1 * J for J in range(1, 9)], [1 / J for J in range(1, 9)], 3.0), gap,
           tone([f2 * J for J in range(1, 13)], [1 / J for J in range(1, 13)], 3.0), gap,
           tone(obs + [third], [1, 1, 1], 3.0)]
    y = np.concatenate(seq)
    y = 0.85 * y / np.abs(y).max()
    wavfile.write(path, sr, (y * 32767).astype(np.int16))
    return dict(observed_Hz=[round(x, 1) for x in obs], third_Hz=round(third, 1), k1_f0_Hz=round(float(f1), 1),
                k2_f0_Hz=round(float(f2), 1), seconds=round(len(y) / sr, 1))


if __name__ == '__main__':
    (a, b), fam = family(*OBS)
    print(f'HerBS-108 lines {OBS} GHz: ratio {OBS[1] / OBS[0]:.5f} ~ {b}/{a}')
    for row in fam:
        print(f"  k={row['k']}: CO J = {row['J']}, z = {row['z']:.4f}, fundamental = {row['f0_GHz']:.3f} GHz")
    print('  (z-GAL quotes z = 1.25 for the CO(2-1)/(3-2) reading; the frequencies imply the k=1 value above)')
    print('\nAnharmonic tie-breaker (two lines read at the wrong octave):')
    tb = tie_breaker()
    for t in tb:
        print(f"  {t['pair']} vs {t['alt']}: {t['mismatch_kms']:.1f} km/s")
    sig30 = np.sqrt(2) * 450 / 2.3548 / 30
    print(f'  CRB centroid precision for FWHM 450 km/s: sqrt(2)*sigma/SNR = {sig30:.1f} km/s per line at S/N 30 '
          f'({np.sqrt(2) * sig30:.1f} km/s on a pair difference)')
    plot(fam, FIG / 'redshift_ambiguity.png')
    info = sonify(str(AUDIO / 'herbs108_octave_ambiguity.wav'))
    print('\naudio:', info)
    (RESULTS / 'redshift_ambiguity.json').write_text(json.dumps(dict(family=fam, tie_breaker=tb, audio=info), indent=1))
