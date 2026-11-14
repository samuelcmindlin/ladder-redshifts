"""Regression tests for the ladder-redshift engine, on synthetic and real ALMA spectra.

Run from the repository root:  python -m pytest   or   python tests/test_ladder_redshifts.py
"""
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ladder_redshifts.ladder import CO, Grid, Engine, sled_flux, simulate, matched_filter, templates_grid, TOL_LN, line_list  # noqa: E402

T_NAMES = line_list()[0]


def test_co_rest_frequencies_match_rotational_constants():
    B, D, H = 57635.968, 0.183505, 1.7e-7            # MHz (ground-state rotational constants of 12CO)
    for J in range(1, 14):
        nu = (2 * B * J - 4 * D * J ** 3 + H * J ** 3 * ((J + 1) ** 3 - (J - 1) ** 3)) / 1e3
        assert abs(nu - CO[J - 1]) < 2e-5, J               # agree to better than 20 kHz


def test_octave_family_of_a_3_to_2_pair():
    from ladder_redshifts.ambiguity import family
    (a, b), fam = family(97.75, 146.60)
    assert (a, b) == (2, 3)
    zs = [r['z'] for r in fam]
    assert abs(zs[1] - 3.7166) < 0.002 and abs(zs[0] - 1.3586) < 0.002
    assert abs((1 + zs[1]) / (1 + zs[0]) - 2) < 1e-3          # candidate redshifts: an overtone series in 1+z


def test_engine_recovers_noiseless_three_line_source():
    class Zero:
        def normal(self, mu, s, n):
            return np.zeros(n)
    g = Grid(66, 164); cov = g.mask([(84.2, 114.9), (125.0, 163.0)])
    f = sled_flux(3.7, 0.29); amp = 10 * f / f[3]
    T, lw = templates_grid()
    e = Engine(g, cov, T, logw=lw)
    Y = matched_filter(simulate(g, 3.0, amp, 450, Zero()), cover=cov)
    z, c, _ = e.posterior(e.loglike(Y))
    assert abs(np.log1p(z) - np.log1p(3.0)) < TOL_LN and c > 0.99


def test_known_traps_flagged_then_resolved():
    from ladder_redshifts import traps
    res, _ = traps.run(nreal=20, ROBUST=True)
    by = {r['name']: r for r in res}
    for r in res:
        for v in r['variants']:
            last = v['steps'][-1]
            assert last['confident_correct'] >= 0.9, (r['name'], v['variant'])     # resolved by the decisive data
    # before the decisive data: never confidently wrong in most draws, and single-line SPT0551-50 stays ambiguous
    for v in by['SPT0551-50']['variants']:
        assert v['steps'][0]['flagged_ambiguous'] >= 0.9
    for r in res:
        for v in r['variants']:
            assert v['steps'][0]['confident_wrong'] <= 0.15, (r['name'], v['variant'])


def test_real_spt0027_band2_spectrum():
    from ladder_redshifts import realspec
    from ladder_redshifts.batch_real import read_spectrum, to_grid_auto
    from ladder_redshifts.ladder import prior_spt
    spec = read_spectrum(ROOT / 'data' / 'SPT0027-50_band2_spectrum.csv')
    r = realspec.run_engine(*to_grid_auto(spec), prior=prior_spt)
    z, c = r['sled']
    assert abs(np.log1p(z) - np.log1p(3.444)) < TOL_LN and c >= 0.95


def test_real_spt0314_admit_single_line_is_ambiguous():
    """ADMIT-route spectrum of a single-line source: the engine must not be confident, must favour CO(3-2)
    (Reuter et al. 2020's identification, since confirmed by Riechers 2025), and must name CO(2-1) at z = 1.62
    as the alternative.
    Without the zero-filled non-detections the candidates spread over at least four rungs."""
    from ladder_redshifts import realspec
    from ladder_redshifts.batch_real import read_spectrum, to_grid_auto
    from ladder_redshifts.ladder import prior_spt
    spec = read_spectrum(ROOT / 'data' / 'admit' / 'spectra' / 'SPT0314-44.csv')
    r = realspec.run_engine(*to_grid_auto(spec), prior=prior_spt)
    z, c = r['sled']
    assert abs(np.log1p(z) - np.log1p(2.933)) < TOL_LN and 0.6 < c < 0.95, (z, c)
    zs = [zz for zz, m in r['candidates'] if m > 0.05]
    assert any(abs(np.log1p(zz) - np.log1p(1.622)) < TOL_LN for zz in zs), r['candidates']
    meas = [w for w in spec if w['kind'] == 'measured']
    rm = realspec.run_engine(*to_grid_auto(meas), prior=prior_spt)
    assert sum(m > 0.05 for _, m in rm['candidates']) >= 4 and rm['sled'][1] < 0.5, rm['candidates']


def test_broad_line_not_split_into_blend_pair():
    """SPT0054-41 (ADMIT route) has one flat-topped ~1100 km/s line, CO(5-4) at z = 4.877. Before the filter-overlap
    term in the likelihood normalization, the engine read it as CO(7-6) + [CI](2-1) at z = 7.23 with confidence
    1.00, counting the flux between the two predicted positions twice."""
    from ladder_redshifts import realspec
    from ladder_redshifts.batch_real import read_spectrum, to_grid_auto
    from ladder_redshifts.ladder import prior_spt, CONF_WIN
    spec = read_spectrum(ROOT / 'data' / 'admit' / 'spectra' / 'SPT0054-41.csv')
    for sub in (spec, [w for w in spec if w['kind'] == 'measured']):
        r = realspec.run_engine(*to_grid_auto(sub), prior=prior_spt)
        p_blend = r['post'][np.abs(r['lnz'] - np.log1p(7.23)) < 0.01].sum()
        assert r['sled'][1] < 0.95 and p_blend < 0.05, (r['sled'], p_blend)
    r = realspec.run_engine(*to_grid_auto(spec), prior=prior_spt)
    assert abs(np.log1p(r['sled'][0]) - np.log1p(4.8769)) < TOL_LN, r['candidates']


def _score_dir(sub):
    from ladder_redshifts import realspec
    from ladder_redshifts.batch_real import read_spectrum, to_grid_auto, read_truth, PRIORS
    truth = read_truth()
    out = {}
    for f in sorted((ROOT / 'data' / sub).glob('*.csv')):
        t = truth.get(f.stem)
        spec = read_spectrum(f)
        if not t or not any(w['kind'] == 'measured' for w in spec):
            continue
        z, c = realspec.run_engine(*to_grid_auto(spec), prior=PRIORS[t['prior']])['sled']
        out[f.stem] = (abs(np.log1p(z) - np.log1p(float(t['z_published']))) < TOL_LN, c, z, t['notes'])
    return out


def test_real_batch_has_no_confident_errors():
    """Every real source with a published z and an on-source spectrum, cube and ADMIT routes: never confident
    and wrong."""
    for sub, nmin in (('spectra', 19), ('admit/spectra', 12)):
        res = _score_dir(sub)
        assert len(res) >= nmin, (sub, len(res))
        for src, (ok, c, z, _) in res.items():
            assert ok or c < 0.95, (sub, src, z, c)


def test_real_cubes_multi_line_sources_confident_single_line_not():
    """Pipeline-cube cut-outs of 2016.1.00672.S: the 7 sources with two or three lines in the published 3 mm scan
    are confident and correct; the 11 with a single line are never confident."""
    res = _score_dir('spectra')
    multi = {s: v for s, v in res.items() if '(2 lines)' in v[3] or '(3 lines)' in v[3]}
    single = {s: v for s, v in res.items() if '(1 line)' in v[3] or ' only;' in v[3]}
    assert len(multi) == 7 and all(ok and c >= 0.95 for ok, c, _, _ in multi.values()), multi
    assert len(single) == 11 and all(c < 0.95 for _, c, _, _ in single.values()), single


# ---- invariances and robustness (added after the pre-publication review, docs/review-response.md) ----

def test_gridding_is_independent_of_window_order():
    from ladder_redshifts import realspec
    from ladder_redshifts.batch_real import read_spectrum, to_grid_auto
    from ladder_redshifts.ladder import prior_spt
    rng = np.random.default_rng(0)
    for src in ('SPT2349-52', 'SPT2152-40', 'SPT2311-45'):
        spec = read_spectrum(ROOT / 'data' / 'spectra' / f'{src}.csv')
        ref = to_grid_auto(spec)
        r0 = realspec.run_engine(*ref, prior=prior_spt)['sled']
        for order in (spec[::-1], [spec[i] for i in rng.permutation(len(spec))]):
            g = to_grid_auto(order)
            assert all(np.array_equal(a, b) for a, b in zip(ref[1:], g[1:])), src
            assert realspec.run_engine(*g, prior=prior_spt)['sled'] == r0, src


def test_overlapping_windows_are_combined_by_inverse_variance():
    from ladder_redshifts.realspec import grid_windows
    f = 100.0 * np.exp(np.arange(40) * 50.0 / 299792.458)
    a = dict(spw='a', freq=f, S=np.full(40, 2.0), sig=np.full(40, 1.0), cont=np.zeros(40))
    b = dict(spw='b', freq=f, S=np.full(40, 4.0), sig=np.full(40, 2.0), cont=np.zeros(40))
    _, d, sig, cov = grid_windows([a, b], Grid(99.0, 102.0))
    s_comb = (2.0 / 1 + 4.0 / 4) / (1 / 1 + 1 / 4)                 # inverse-variance mean flux
    n_comb = (1 / 1 + 1 / 4) ** -0.5
    assert np.allclose(sig[cov], n_comb) and np.allclose(d[cov], s_comb / n_comb)


def test_gaps_inside_a_window_are_not_bridged():
    from ladder_redshifts.realspec import grid_windows
    f = np.r_[100.0 * np.exp(np.arange(10) * 50.0 / 299792.458), 101.0 * np.exp(np.arange(10) * 50.0 / 299792.458)]
    w = dict(spw='w', freq=f, S=np.zeros(20), sig=np.ones(20), cont=np.zeros(20))
    grid, _, _, cov = grid_windows([w], Grid(99.0, 102.0))
    assert cov.sum() <= 22                                           # ~10 pixels per fragment, not ~200


def test_data_outside_the_coverage_have_no_effect():
    grid = Grid(90, 110)
    cover = grid.mask([(95, 100)])
    T, lw = templates_grid()
    e = Engine(grid, cover, T, logw=lw)
    rng = np.random.default_rng(3)
    d1 = rng.normal(0, 1, grid.npix)
    d2 = d1.copy()
    d2[~cover] = rng.normal(0, 100, (~cover).sum())                  # change every unobserved pixel
    d2[np.flatnonzero(cover)[-1] + 1] = 1e4
    assert np.array_equal(e.loglike(e.filter(d1)), e.loglike(e.filter(d2)))


def test_filter_output_has_unit_noise_at_window_edges():
    grid = Grid(90, 110)
    cover = grid.mask([(95, 100)])
    edge = np.flatnonzero(cover)[0]
    rng = np.random.default_rng(4)
    y = np.array([matched_filter(rng.normal(0, 1, grid.npix), cover=cover)[:, [edge, edge + 3]] for _ in range(600)])
    assert np.all(np.abs(y.std(axis=0) - 1) < 0.1), y.std(axis=0)


def test_equal_weight_baseline_weights_line_snr_equally():
    """The comparison baseline must not inherit the noise weighting of the SLED likelihood."""
    from ladder_redshifts import realspec
    from ladder_redshifts.batch_real import read_spectrum, to_grid_auto
    spec = read_spectrum(ROOT / 'data' / 'spectra' / 'SPT0436-40.csv')
    spec = [dict(w, sig=w['sig'] * (5.0 if i % 2 else 1.0), S=w['S'] * (5.0 if i % 2 else 1.0))
            for i, w in enumerate(spec)]                             # very unequal noise between windows
    grid, d, sig, cov = to_grid_auto(spec)
    r = realspec.run_engine(grid, d, sig, cov)
    e = Engine(grid, cov, np.ones(len(T_NAMES)), zmin=0.5, zmax=8.0)
    Y = e.filter(d)
    z, c, _ = e.posterior(e.loglike(Y), refine=e.loglike(Y, dv_max=0))
    assert r['flat'] == (float(z), float(c))


def test_detection_width_is_chosen_by_value():
    from ladder_redshifts import realspec
    from ladder_redshifts.batch_real import read_spectrum, to_grid_auto
    from ladder_redshifts.ladder import BANK
    g = to_grid_auto(read_spectrum(ROOT / 'data' / 'spectra' / 'SPT0403-58.csv'))
    ref = realspec.run_engine(*g)
    one = realspec.run_engine(*g, bank=(450.0,))                     # used to raise IndexError
    rev = realspec.run_engine(*g, bank=tuple(reversed(BANK)))
    assert one['detections'] == ref['detections'] == rev['detections'] and rev['sled'][0] == ref['sled'][0]


def test_extraction_failure_gives_nonzero_exit_status():
    try:
        from astropy.io import fits
    except ImportError:                                              # astropy is only needed for scripts/alma
        return
    import subprocess
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        cube = Path(tmp) / 'x.TEST_sci.spw0.cube.I.pbcor.fits'
        fits.PrimaryHDU().writeto(cube)                              # no frequency axis
        for extra in ([], ['--tag', 'freq']):
            proc = subprocess.run([sys.executable, str(ROOT / 'scripts' / 'alma' / 'extract_cube_spectra.py'), str(cube),
                                   '--out', str(Path(tmp) / 'out'), '--scratch', str(Path(tmp) / 's')] + extra,
                                  capture_output=True, text=True)
            assert proc.returncode != 0 and 'FAILED' in proc.stdout, (extra, proc.returncode, proc.stdout[-300:])


# ---- noise scale, velocity offsets, shared harmonics (added 12 November 2026; docs/development-notes.md, 20-23) ----

def _cube_route_files():
    return sorted((ROOT / 'data' / 'spectra').glob('*.csv')) + sorted((ROOT / 'data').glob('*_band2_spectrum*.csv'))


def test_shipped_cube_spectra_have_unit_noise_scatter():
    """The line-free channels of flux/noise must scatter by 1. The extractor's annulus noise was about 1.3 times
    too high at the source (1.7 for SPT0027-50), which hid 5-7 sigma lines from the detection rule."""
    from ladder_redshifts.batch_real import read_spectrum
    from ladder_redshifts.realspec import noise_scale
    files = _cube_route_files()
    assert len(files) >= 24
    for f in files:
        sc = np.array([noise_scale(w['S'], w['sig']) for w in read_spectrum(f)])
        assert abs(np.median(sc) - 1) < 0.05, (f.name, np.median(sc))
        assert sc.min() > 0.7 and sc.max() < 1.45, (f.name, sc.min(), sc.max())


def test_calibrate_noise_recovers_a_wrong_noise_scale():
    from ladder_redshifts.realspec import calibrate_noise
    rng = np.random.default_rng(6)
    over = 1.25 + 0.03 * np.arange(12)                               # stated noise / true noise, rising with frequency
    spec = []
    for k in range(12):
        f = (86 + 2.4 * k) * np.exp(np.arange(120) * 50.0 / 299792.458)
        S = rng.normal(0, 1.0, 120)
        if k == 4:
            S[50:60] += 25.0                                         # a bright line must not bias the scale
        spec.append(dict(spw=f'w{k}', freq=f, S=S, sig=np.full(120, over[k]), cont=np.zeros(120)))
    cal = calibrate_noise(spec)
    applied = np.array([w['noise_scale'] for w in cal])
    assert np.all(np.abs(applied * over - 1) < 0.08), applied * over
    assert all(np.allclose(w['sig'], 1.0, rtol=0.08) for w in cal)
    zero = dict(spw='z', freq=spec[0]['freq'], S=np.zeros(120), sig=np.full(120, 10.0), cont=np.zeros(120),
                kind='zero-filled')                                  # ADMIT non-detections are not spectra
    assert np.array_equal(calibrate_noise([zero])[0]['sig'], zero['sig'])


def test_point_source_spectrum_has_calibrated_noise():
    """spectrum_point() on synthetic cut-outs whose annulus noise is 1.4 times the noise at the source."""
    from ladder_redshifts.realspec import spectrum_point, noise_scale
    rng = np.random.default_rng(8)
    spws = []
    for k in range(8):
        freq = (88 + 3.0 * k) * np.exp(np.arange(119) * 50.0 / 299792.458)
        cut = rng.normal(0, 1e-3, (119, 1, 1)) * np.ones((1, 9, 9))  # Jy/beam; one beam covers the cut-out
        spws.append(dict(name=f'T_f{k}', freq=freq, cut=cut, rms=np.full(119, 1.4e-3),
                         beams=np.array([[4.0, 3.5, 0.0]]), pix=0.45, k=2))
    spec = spectrum_point(spws)
    assert abs(np.median([w['noise_scale'] for w in spec]) * 1.4 - 1) < 0.07
    assert abs(np.median([noise_scale(w['S'], w['sig']) for w in spec]) - 1) < 0.1


class _NoNoise:
    def normal(self, mu, s, n):
        return np.zeros(n)


def test_line_offsets_within_the_tolerance_leave_the_likelihood_unchanged():
    """Two lines, one displaced by 100 km/s. With the per-line tolerance the likelihood of the solution does not
    change; with exact positions (dv_max = 0) the same offset costs likelihood, which at high S/N is enough to
    decide between readings that the offset should not distinguish."""
    g = Grid(66, 164); cov = g.mask([(84.2, 114.9), (125.0, 163.0)])
    T, lw = templates_grid()
    e = Engine(g, cov, T, logw=lw)
    f = sled_flux(3.7, 0.29); amp = 30 * f / f[3]
    j = T_NAMES.index('CO(5-4)')
    res = {}
    for dv in (-100, 0, 100):
        rest, moved = amp.copy(), np.zeros(len(amp))
        moved[j], rest[j] = amp[j], 0.0
        d = simulate(g, 3.2, rest, 450, _NoNoise()) + simulate(g, 4.2 * np.exp(dv / 299792.458) - 1, moved, 450, _NoNoise())
        Y = e.filter(d)
        L, L0 = e.loglike(Y), e.loglike(Y, dv_max=0)
        z, c, _ = e.posterior(L, refine=L0)
        assert abs(np.log1p(z) - np.log1p(3.2)) < TOL_LN and c > 0.99, (dv, z, c)
        res[dv] = (L.max(), L0.max())
    for dv in (-100, 100):
        assert abs(res[dv][0] / res[0][0] - 1) < 0.02, res         # tolerant: unchanged
        assert res[dv][1] < 0.985 * res[0][1], res                  # exact positions: penalised


def test_line_predicted_in_a_one_pixel_gap_still_counts():
    """A reading must not escape the penalty for a line it predicts but the data lack by parking that line in
    the one-pixel gap between two spectral windows (it can shift by a pixel or two within the tolerance)."""
    g = Grid(66, 164)
    cov = g.mask([(84.2, 114.9)])
    T, lw = templates_grid()
    f = sled_flux(3.7, 0.29)
    d = simulate(g, CO[2] / 90.0 - 1, np.where(np.arange(len(f)) == 2, 30.0, 0.0), 450, _NoNoise())   # one line, 90 GHz
    z_alt = CO[4] / 90.0 - 1                                         # read as CO(5-4): CO(6-5) should be at 108 GHz
    mass = []
    for hole in (False, True):
        c = cov.copy()
        if hole:
            p = int(round(g.pix(CO[5] / (1 + z_alt))))
            c[[p - 2, p + 2]] = False                                # gaps within reach of the predicted pixel
        e = Engine(g, c, T, logw=lw)
        _, _, post = e.posterior(e.loglike(e.filter(d)))
        mass.append(post[np.abs(e.lnz - np.log1p(z_alt)) < 0.006].sum())
    assert mass[0] < 1e-6 and mass[1] < 1e-6, mass


def test_shared_harmonic_answer_does_not_depend_on_velocity_offsets():
    """SPT0314-44 plus a simulated tuning that covers only the harmonic both candidates share (CO(9-8) at
    z = 2.933 is CO(6-5) at z = 1.622, 76 km/s apart through the ladder's anharmonicity). Before the tolerance,
    displacing that line by 25-76 km/s turned a correct answer at confidence 1.00 into a wrong one at 1.00.
    A tuning with lines that only one candidate predicts must still settle the redshift at any offset."""
    from ladder_redshifts import followup as F
    from ladder_redshifts import realspec
    from ladder_redshifts.batch_real import read_spectrum, to_grid_auto
    from ladder_redshifts.ladder import prior_spt
    spec = read_spectrum(ROOT / 'data' / 'spectra' / 'SPT0314-44.csv')
    grid, d, sig, cov = to_grid_auto(spec)
    r0 = realspec.run_engine(grid, d, sig, cov, prior=prior_spt)
    f_det, snr = max(r0['detections'], key=lambda q: q[1])
    amp = snr * float(sig[int(round(grid.pix(f_det)))])              # flux of the observed line
    sig_new = float(np.median(np.concatenate([w['sig'] for w in spec])))
    covered = sorted((float(w['freq'].min()), float(w['freq'].max())) for w in spec)

    def answer(lo, band, z_true, line, dv):
        wins = F.tuning(lo, band, covered)
        z_new = (1 + z_true) * np.exp(-dv / 299792.458) - 1          # the new lines, displaced by dv km/s
        new = F.synth_windows(wins, z_new, amp / sig_new, T_NAMES.index(line), 3.7, sig_new, 400.0,
                              np.random.default_rng(1))
        return realspec.run_engine(*to_grid_auto(spec + new), prior=prior_spt)['sled']

    for z_true, line in ((2.9332, 'CO(3-2)'), (1.6222, 'CO(2-1)')):
        z0, c0 = answer(255.0, 'Band 6', z_true, line, 0)            # shared harmonic only
        for dv in (-100, 100):
            z, c = answer(255.0, 'Band 6', z_true, line, dv)
            assert abs(np.log1p(z) - np.log1p(z0)) < TOL_LN and abs(c - c0) < 0.04, (z_true, dv, (z, c), (z0, c0))
        for dv in (-100, 100):                                       # Band 4: CO(5-4) or CO(3-2), one per candidate
            z, c = answer(139.0, 'Band 4', z_true, line, dv)
            assert abs(np.log1p(z) - np.log1p(z_true)) < TOL_LN and c >= 0.95, (z_true, dv, z, c)


def test_planner_does_not_count_tunings_that_only_add_shared_lines():
    from ladder_redshifts import followup as F
    nus = line_list()[1]
    covered = [(84.1, 114.8)]
    assert not F.separates(F.tuning(255.0, 'Band 6', covered), 2.9332, 1.6222, nus)            # CO(9-8) = CO(6-5)
    assert F.separates(F.tuning(139.0, 'Band 4', covered), 2.9332, 1.6222, nus)                # CO(5-4) vs CO(3-2)


def test_follow_up_tunings_fit_the_alma_receivers():
    """Every tuning the planner offers is two 3.75 GHz sidebands inside the band, at an offset from the LO that
    the band's IF range allows (4-8 GHz; 4.5-10 GHz in Band 6). The earlier planner put them at 6-9.75 GHz in
    every band, which no 4-8 GHz receiver can do."""
    from ladder_redshifts import followup as F
    for b, (f_lo, f_hi) in F.BANDS.items():
        near, far = F.sideband_offsets(b)
        assert F.IF_RANGE[b][0] <= near < far <= F.IF_RANGE[b][1] and abs(far - near - F.SIDEBAND_GHZ) < 1e-9
        los = F.lo_grid(b)
        assert len(los) > 10 and los[0] - far >= f_lo and los[-1] + far <= f_hi
        assert los[0] - 1.0 - far < f_lo and los[-1] + 1.0 + far > f_hi              # the grid spans the whole range
        for lo in los:
            wins = F.tuning(lo, b, [])
            assert len(wins) == 2 and all(abs(y - x - F.SIDEBAND_GHZ) < 1e-9 for x, y in wins)
            assert f_lo <= wins[0][0] and wins[1][1] <= f_hi
    assert F.sideband_offsets('Band 3') == (4.125, 7.875) and F.sideband_offsets('Band 6') == (5.375, 9.125)


# ---- comparisons added on 12 November 2026: guided rule, degradation test, the Sudoku code's classes ----

def test_guided_rule_confirms_a_weak_line_and_otherwise_equals_the_detection_rule():
    """CO(4-3) at S/N 20 with [CI](1-0) at S/N 4 (z = 3.85, Band 3 scan), no noise. The detection rule sees one
    line and cannot answer; the guided rule finds the second line at the position z = 3.85 predicts, provided its
    threshold is below that line. With both lines above 5 sigma the two rules agree."""
    from ladder_redshifts.ladder import detect, sudoku_rule, guided_rule, detection_index, BANK

    class Zero:
        def normal(self, mu, s, n):
            return np.zeros(n)
    g = Grid(80, 119); cov = g.mask([(84.2, 114.9)])
    names, nus = line_list()
    amp = np.zeros(len(nus)); amp[names.index('CO(4-3)')] = 20.0

    def rules(ci):
        amp[names.index('[CI](1-0)')] = ci
        y = matched_filter(simulate(g, 3.85, amp, 450, Zero()), cover=cov)[detection_index(BANK)]
        dets = detect(y, cov)
        return dets, sudoku_rule(dets, g, cov, y)[0], [guided_rule(dets, g, cov, y, confirm=t)[0] for t in (3.0, 4.5)]

    dets, z_rule, (z_g3, z_g45) = rules(4.0)
    assert len(dets) == 1 and z_rule is None and z_g45 is None
    assert z_g3 is not None and abs(np.log1p(z_g3) - np.log1p(3.85)) < TOL_LN
    dets, z_rule, (z_g3, z_g45) = rules(8.0)
    assert len(dets) == 2 and z_rule is not None and z_rule == z_g3 == z_g45


def test_degradation_reaches_the_requested_second_line_snr():
    """Adding noise to SPT0652-55 until its [CI](1-0) line should sit at S/N 4: the measured S/N at the line is
    about 4, and every method's outcomes add up."""
    from ladder_redshifts import degrade as D
    from ladder_redshifts.batch_real import read_spectrum, to_grid_auto
    from ladder_redshifts.realspec import run_engine
    from ladder_redshifts.ladder import prior_spt
    src, z = 'SPT0652-55', 3.3466
    assert (src, z) in D.multi_line_sources() and len(D.multi_line_sources()) == 7
    r = run_engine(*to_grid_auto(read_spectrum(ROOT / 'data' / 'spectra' / f'{src}.csv')), prior=prior_spt)
    lines = D.published_lines(r['detections'], z)
    assert [n for n, _, _ in lines] == ['CO(4-3)', '[CI](1-0)']
    _, f2, snr2 = lines[1]
    _, out = D.run_level((src, z, f2, 4.0, snr2 / 4.0, 12, [1, 2, 3]))
    assert abs(out['second_line_snr_measured'] - 4.0) < 0.8
    assert abs(sum(out['rule'].values()) - 1) < 1e-9 and all(abs(sum(v.values()) - 1) < 1e-9 for v in out['guided'].values())
    e = out['spt']['engine']
    assert e['confident_correct'] + e['confidently_wrong'] <= 1 and e['best_is_published'] >= e['confident_correct']
    # the source's own photometric redshift, in its three versions, is run as a prior as well
    from ladder_redshifts.photoz import STATES
    assert set(out['guided_photoz']) == set(STATES) <= set(out) and set(STATES) <= set(out['xcorr'])
    assert all(abs(sum(v.values()) - 1) < 1e-9 for v in out['guided_photoz'].values())
    assert set(out['guided_co_only']) == set(STATES) | {'spt', 'uniform'}
    assert all(abs(sum(v.values()) - 1) < 1e-9 for v in out['guided_co_only'].values())
    # every setting of the rule family, and the engine at each confidence threshold, are counted per draw
    fam, curve = out['rule_family'], out['engine_curve']
    assert [len(fam[p]) for p in D.FAMILY_PRIORS] == [72, 216] and all(0 <= a + b <= 12 for a, b in fam['photoz'])
    j = D.rule_family().index(dict(detect=5.0, confirm=3.0, ci_counted=True, prune_expect=8.0, share=None))
    assert fam['uniform'][j] == [round(out['guided']['3']['correct'] * 12), round(out['guided']['3']['wrong'] * 12)]
    k = D.ENGINE_THRESHOLDS.index(0.95)
    assert curve['uniform'][k] == [round(out['uniform']['engine']['confident_correct'] * 12),
                                   round(out['uniform']['engine']['confidently_wrong'] * 12)]
    assert all(curve['uniform'][i][0] >= curve['uniform'][i + 1][0] for i in range(len(D.ENGINE_THRESHOLDS) - 1))


def test_sudoku_comparison_file_matches_the_current_engine_results():
    """results/redshift_sudoku_comparison.json quotes this engine's verdicts next to the classes returned by the
    High-z Sudoku code (scripts/compare_sudoku.py, which needs that code and is not run here). The quoted
    verdicts must be the current ones, and the classes must cover every source with a spectrum and a redshift."""
    import csv
    import json
    o = json.loads((ROOT / 'results' / 'redshift_sudoku_comparison.json').read_text())
    batch = {r['source']: r for r in json.loads((ROOT / 'results' / 'redshift_real_batch.json').read_text())}
    with open(ROOT / 'data' / 'truth.csv', newline='') as fh:
        want = {r['source'] for r in csv.DictReader(fh) if (ROOT / 'data' / 'spectra' / f"{r['source']}.csv").exists()}
    assert {r['source'] for r in o['sources']} == want and len(want) == 19
    for r in o['sources']:
        b = batch[r['source']]
        assert r['engine'] == dict(z=b['z_engine'], confidence=b['confidence'], verdict=b['verdict']), r['source']
        assert r['sudoku_co']['lines'] <= r['sudoku_co_ci']['lines'] <= 3


def test_guided_rule_uses_a_prior_only_when_one_solution_dominates():
    """One line, CO(3-2) at S/N 20 (z = 2.64), no noise: the line leaves two solutions (CO(2-1) at z = 1.43 and
    CO(3-2) at z = 2.64) and nothing confirms either. Without a prior the guided rule cannot answer. The broad
    survey prior does not single one out either (0.91 of the weight, below the 0.95 required). A prior that puts
    nearly all its weight on one solution does."""
    from ladder_redshifts.ladder import detect, guided_rule, detection_index, prior_spt, BANK

    class Zero:
        def normal(self, mu, s, n):
            return np.zeros(n)
    g = Grid(80, 119); cov = g.mask([(84.2, 114.9)])
    names, nus = line_list()
    amp = np.zeros(len(nus)); amp[names.index('CO(3-2)')] = 20.0
    y = matched_filter(simulate(g, 2.64, amp, 450, Zero()), cover=cov)[detection_index(BANK)]
    dets = detect(y, cov)
    z0, nsol, nconf = guided_rule(dets, g, cov, y)
    assert z0 is None and nsol == 2 and nconf == 0
    assert guided_rule(dets, g, cov, y, prior=prior_spt)[0] is None
    z1 = guided_rule(dets, g, cov, y, prior=lambda z: np.exp(-0.5 * ((z - 2.64) / 0.05) ** 2))[0]
    assert z1 is not None and abs(np.log1p(z1) - np.log1p(2.64)) < TOL_LN


def test_cross_correlation_noise_matches_shuffles_and_two_lines_give_the_redshift():
    """The template cross-correlation (after Yun et al. 2015). (1) Its closed-form noise is the limit of the
    paper's shuffling: with 1,500 shuffles the two S/N values agree to a few per cent wherever the S/N is
    large. (2) With CO(4-3) at S/N 15 and [CI](1-0) at S/N 6 in a Band 3 scan it returns z = 3.85. (3) With
    CO(4-3) alone the readings of that line as CO(3-2) and as CO(4-3) score the same to a few per cent, so the
    statistic cannot choose between them; multiplied by the survey prior, which peaks at z = 3.9, it does."""
    from ladder_redshifts.ladder import cross_correlation, shuffled_cross_correlation, xcorr_redshift, prior_spt
    g = Grid(80, 119); cov = g.mask([(84.2, 114.9)])
    names, nus = line_list()
    amp = np.zeros(len(nus)); amp[names.index('CO(4-3)')] = 15.0; amp[names.index('[CI](1-0)')] = 6.0
    d = simulate(g, 3.85, amp, 450, np.random.default_rng(5))
    lnz, snr = cross_correlation(d, g, cov)
    _, snr_s = shuffled_cross_correlation(d, g, cov, 1500, np.random.default_rng(6))
    big = np.abs(snr) > 4
    assert big.sum() > 20 and np.all(np.abs(snr_s[big] / snr[big] - 1) < 0.12)
    assert abs(np.median(snr_s[big] / snr[big]) - 1) < 0.04
    z, s_best = xcorr_redshift(d, g, cov)
    assert abs(np.log1p(z) - np.log1p(3.85)) < TOL_LN and s_best > 10
    amp[names.index('[CI](1-0)')] = 0.0
    d1 = simulate(g, 3.85, amp, 450, np.random.default_rng(5))
    lnz, snr1 = cross_correlation(d1, g, cov)
    one_line = [np.max(snr1[np.abs(lnz - (np.log(nus[j]) - np.log(nus[3]) + np.log(4.85))) < TOL_LN]) for j in (2, 3)]
    assert abs(one_line[0] - one_line[1]) < 0.06 * one_line[1]           # CO(3-2) at z = 2.64 against CO(4-3) at 3.85
    zp, _ = xcorr_redshift(d1, g, cov, prior=prior_spt)
    assert abs(np.log1p(zp) - np.log1p(3.85)) < TOL_LN                   # the survey prior peaks at z = 3.9


if __name__ == '__main__':
    import inspect
    fails = 0
    for name, fn in sorted(globals().items()):
        if name.startswith('test_') and inspect.isfunction(fn):
            try:
                fn(); print('PASS', name)
            except AssertionError as e:
                fails += 1; print('FAIL', name, e)
    sys.exit(1 if fails else 0)


def test_photometric_redshift_table_matches_its_calibration_record():
    """data/photoz.csv lists every SPT source of the truth table, and the recalibration used by the code is the one
    scripts/photoz_calibration.py derived from the 60 SPT sources outside this sample."""
    import json
    from ladder_redshifts import photoz as P
    from ladder_redshifts.batch_real import read_truth
    table = P.read_photoz()
    cal = json.loads((ROOT / 'results' / 'redshift_photoz_calibration.json').read_text())
    assert sorted(table) == cal['sample'] == sorted(s for s in read_truth() if s.startswith('SPT')) and len(table) == 21
    assert all(0.2 <= e <= 1.5 and 0.2 <= e2 <= 1.5 and 1.5 < z < 7 and 1.5 < z2 < 7 for z, e, z2, e2 in table.values())
    rc = cal['recalibration']
    assert cal['published_2020']['outside_sample']['n'] == 60
    assert (P.RECAL['shift'], P.RECAL['width']) == (rc['shift_in_published_errors'], rc['width_factor'])
    st = P.photoz_states('SPT0314-44')                       # 2.8 +/- 0.5 in the 2020 paper, 2.9 +/- 0.4 in 2024
    assert st['photoz'] == (2.8, 0.5) and st['photoz_2024'] == (2.9, 0.4)
    assert abs(st['photoz_recal'][0] - (2.8 - 0.479 * 0.5)) < 1e-12 and abs(st['photoz_recal'][1] - 0.5505) < 1e-12
    assert P.photoz_states('HerBS-108') == {}
    g, w = P.gaussian(3.0, 0.5), P.window(3.0)
    assert abs(g(3.5) - np.exp(-0.5)) < 1e-12 and w(3.5) == 1 and w(3.53) == 0


def test_photometric_redshift_priors_settle_single_line_sources_without_a_wrong_call():
    """With each source's published photometric redshift as the prior, the engine becomes confident of some of the
    one-line sources and is never confident of a wrong redshift; the saved results are those of the current code."""
    import json
    from ladder_redshifts import photoz as P
    rows = P.analyse()
    saved = json.loads((ROOT / 'results' / 'redshift_photoz.json').read_text())
    assert json.loads(json.dumps(rows)) == saved['sources']
    assert json.loads(json.dumps(P.summarise(rows))) == saved['summary']
    one, several = saved['summary']['one'], saved['summary']['several']
    assert (one['n'], several['n']) == (11, 8)
    assert one['reference']['spt']['confident_correct'] == 0 and one['reference']['spt']['best_is_published'] == 10
    for g in (one, several):
        assert all(s['engine']['confident_wrong'] == 0 for s in g['states'].values())
    assert [one['states'][s]['engine']['confident_correct'] for s in P.STATES] == [7, 6, 4]
    assert all(s['engine']['confident_correct'] == 8 for s in several['states'].values())
    # SPT0528-53: its [CI] line is missing, so the rule drops the published redshift unless only CO lines count
    r = next(r for r in rows if r['source'] == 'SPT0528-53')
    near = lambda zs: any(abs(np.log1p(z) - np.log1p(r['z_published'])) < TOL_LN for z in zs)
    assert near(r['rule_solutions']['co_only']) and not near(r['rule_solutions']['ci_counted'])
    # SPT2349-52: the 2024 photometric redshift is 2.75 sigma from the published redshift and every method follows it
    q = next(r for r in rows if r['source'] == 'SPT2349-52')['states']['photoz_2024']
    assert q['pull'] == 2.75 and q['engine']['outcome'] == 'wrong' and q['engine']['confidence'] < 0.95
    assert q['guided_ci_counted']['outcome'] == q['window_co_only']['outcome'] == q['xcorr']['outcome'] == 'wrong'


def test_the_rule_the_engine_is_compared_with_is_the_best_of_its_family_by_the_stated_criterion():
    """results/redshift_rule_family.csv lists every rule setting; the setting named in the results file is the one
    with the most correct answers among those that are wrong no more often than the engine; the setting used
    before the family existed reproduces the earlier column; and the stored answers on a real spectrum are those
    of the current code."""
    import csv
    import json
    from ladder_redshifts import degrade as D
    from ladder_redshifts.batch_real import read_spectrum, to_grid_auto
    from ladder_redshifts.realspec import engine_likelihood
    from ladder_redshifts.photoz import photoz_priors
    o = json.loads((ROOT / 'results' / 'redshift_degradation.json').read_text())
    with open(ROOT / 'results' / 'redshift_rule_family.csv', newline='') as fh:
        rows = list(csv.DictReader(fh))
    fam = o['rule_family']
    assert fam['draws'] == 5600 and [fam['priors'][p]['n_settings'] for p in D.FAMILY_PRIORS] == [72, 216]
    for p in D.FAMILY_PRIORS:
        rr = [r for r in rows if r['prior'] == p]
        f = fam['priors'][p]
        cap = f['engine_confident_wrong']
        ok = [r for r in rr if int(r['total_wrong']) <= cap]
        best = max(ok, key=lambda r: int(r['total_correct']))
        ch = f['chosen']['no_more_wrong_than_engine']
        assert len(rr) == f['n_settings'] and len(ok) == f['settings_with_no_more_wrong_than_engine']
        assert (ch['total_correct'], ch['total_wrong']) == (int(best['total_correct']), int(best['total_wrong']))
        assert ch['setting'] == D.rule_family(p != 'uniform')[ch['index']] and int(best['index']) == ch['index']
        e95 = next(q for q in f['engine_by_threshold'] if q['threshold'] == 0.95)
        assert e95['engine_wrong'] == cap and e95['rule_index'] == ch['index']
    # the single setting of the earlier comparison is one member of the family
    old = next(r for r in rows if r['prior'] == 'uniform' and (r['detect'], r['confirm'], r['ci_counted'],
                                                               r['prune_expect']) == ('5.0', '3.0', 'True', '8.0'))
    for q in o['summary']:
        assert abs(float(old[f"correct_{q['target']:g}"]) - q['guided']['3']['correct']) < 1e-9
    # with no prior the engine is not ahead of the best rule at its own threshold; with the photometric prior it is
    u, ph = fam['priors']['uniform'], fam['priors']['photoz']
    at = lambda f, t: next(q for q in f['engine_by_threshold'] if q['threshold'] == t)
    assert at(u, 0.95)['rule_correct'] > at(u, 0.95)['engine_correct'] and at(u, 0.8)['engine_correct'] > at(u, 0.8)['rule_correct']
    assert all(q['engine_correct'] > q['rule_correct'] for q in ph['engine_by_threshold'])
    # one real spectrum: the stored letters are the answers of the current code
    src = 'SPT0528-53'
    row = next(r for r in o['as_observed'] if r['source'] == src)
    g = to_grid_auto(read_spectrum(ROOT / 'data' / 'spectra' / f'{src}.csv'))
    like = engine_likelihood(*g)
    for p, prior in (('uniform', None), ('photoz', photoz_priors(src)['photoz'])):
        letters = ''.join('n' if z is None else 'c' if D.correct(z, row['z_published']) else 'w'
                          for z in D.family_answers(like['yd'], g[0], g[3], prior))
        assert letters == row['rule_family'][p]


def test_wrong_rate_bound_and_the_second_sample():
    """The exact bound on the wrong-answer rate; and the independent sample is defined, disjoint from the
    development sample, and reached through $LADDER_REDSHIFTS_SAMPLE without any other change."""
    import csv
    import json
    import os
    import subprocess
    from ladder_redshifts.photoz import upper_bound
    assert abs(upper_bound(0, 15) - (1 - 0.05 ** (1 / 15))) < 1e-9 and abs(upper_bound(0, 8) - 0.312) < 1e-3
    assert upper_bound(0, 59) < 0.05 < upper_bound(0, 58) and upper_bound(1, 15) > upper_bound(0, 15)
    saved = json.loads((ROOT / 'results' / 'redshift_photoz.json').read_text())['summary']['all']['engine']
    assert (saved['spt']['confident_answers'], saved['photoz']['confident_answers']) == (8, 15)
    assert all(v['wrong'] == 0 for v in saved.values()) and saved['photoz']['wrong_rate_upper_95'] == 0.181
    read = lambda rel: list(csv.DictReader(open(ROOT / rel, newline='')))
    new, dev = read('data/independent/truth.csv'), read('data/truth.csv')
    names = [r['source'] for r in new]
    assert len(names) == 23 and not set(names) & {r['source'] for r in dev}
    assert [r['source'] for r in read('data/independent/photoz.csv')] == names
    assert [r['target'] for r in read('data/cutouts/2015.1.00504.S_targets.csv')] == names
    tags = [r['notes'][r['notes'].rindex('('):] for r in new]
    assert (sum(t in ('(2 lines)', '(3 lines)') for t in tags), tags.count('(1 line)')) == (13, 9)
    env = dict(os.environ, LADDER_REDSHIFTS_SAMPLE='independent', LADDER_REDSHIFTS_ROOT=str(ROOT))
    code = ('from ladder_redshifts import paths as p; from ladder_redshifts.batch_real import read_truth; '
            'print(p.DATA.relative_to(p.ROOT).as_posix(), len(read_truth()))')
    out = subprocess.run([sys.executable, '-c', code], env=env, capture_output=True, text=True, cwd=str(ROOT))
    assert out.stdout.split() == ['data/independent', '23'], out.stderr
    # archive names of the re-imaged (ARI-L) cubes that the Cycle 3 project has
    name = 'member.uid___A001_X2f7_X446.ari_l.SPT0020-51_sci.spw0_85009MHz.12m.cube.I.pbcor.fits'
    lst = ROOT / 'tests' / '_ari_l_names.txt'
    man = ROOT / 'tests' / '_ari_l_manifest.tsv'
    try:
        lst.write_text(name + '\n' + name.replace('SPT0020-51', 'SPT0027-50') + '\n')
        proc = subprocess.run([sys.executable, str(ROOT / 'scripts' / 'alma' / 'make_cutout_manifest.py'), '--ari-l',
                               '--files', str(lst), '--targets', str(ROOT / 'data/cutouts/2015.1.00504.S_targets.csv'),
                               '--out', str(man)], capture_output=True, text=True)
        assert proc.returncode == 0, proc.stderr
        body = man.read_text().splitlines()
        assert len(body) == 2 and body[1].split('\t')[:3] == [name, 'cutout', 'SPT0020-51']
    finally:
        for f in (lst, man):
            if f.exists():
                f.unlink()


def test_the_independent_sample_results_are_those_of_the_current_code():
    """The 23 spectra of the independent test, re-run from the saved spectra in a separate process with
    $LADDER_REDSHIFTS_SAMPLE set, give the stored results; and the stored report has the outcomes quoted in the
    README: no confident wrong answer, 9 of 13 multi-line sources settled with the survey prior."""
    import json
    import os
    import subprocess
    code = ('import json; from ladder_redshifts import photoz as P; from ladder_redshifts.paths import RESULTS; '
            'rows = json.loads(json.dumps(P.analyse())); saved = json.loads((RESULTS / "redshift_photoz.json").read_text()); '
            'print(len(rows), rows == saved["sources"], json.loads(json.dumps(P.summarise(rows))) == saved["summary"])')
    env = dict(os.environ, LADDER_REDSHIFTS_SAMPLE='independent', LADDER_REDSHIFTS_ROOT=str(ROOT))
    out = subprocess.run([sys.executable, '-c', code], env=env, capture_output=True, text=True, cwd=str(ROOT))
    assert out.stdout.split() == ['23', 'True', 'True'], out.stderr[-2000:]
    rep = json.loads((ROOT / 'results' / 'independent' / 'report.json').read_text())
    pr, g = rep['primary'], rep['groups']
    assert (pr['spt']['confident_answers'], pr['spt']['wrong'], pr['photoz']['confident_answers'], pr['photoz']['wrong']) \
        == (11, 0, 17, 0)
    assert (pr['spt']['pooled_confident_answers'], pr['photoz']['pooled_confident_answers']) == (19, 32)
    assert pr['photoz']['pooled_wrong_rate_upper_95'] == 0.089
    assert (g['several']['n'], g['one']['n'], g['none']['n']) == (13, 9, 1)
    assert g['several']['engine']['spt']['confident_correct'] == 9 and g['one']['engine']['photoz']['confident_correct'] == 6
    assert g['several']['chosen_rule']['uniform'] == dict(correct=7, wrong=1, none=5)
    d = rep['as_observed_engine_against_chosen_rule']
    assert len(d['photoz']['engine_only']) == 7 and d['photoz']['rule_only'] == [] and d['uniform']['rule_only'] == []
    batch = json.loads((ROOT / 'results' / 'independent' / 'redshift_real_batch.json').read_text())
    assert sum(r['verdict'].startswith('confident & correct') for r in batch) == 11
    assert not any(r['verdict'].startswith('confident') and 'correct' not in r['verdict'] for r in batch)

