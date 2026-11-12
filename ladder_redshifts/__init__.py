"""ladder_redshifts: redshifts from molecular line ladders, treated as missing-fundamental pitch estimation.

Modules
  ladder      the engine: line catalogue, log-frequency grid, simulator, matched filter, estimators
              (SLED-weighted likelihood, equal-weight sum, simplified Sudoku-style rule, naive spacing)
  realspec    real ALMA spectra: aperture / point-source extraction from cube cut-outs; SPT0027-50 (Band 2)
  batch_real  score every real spectrum under data/ against data/truth.csv
  followup    which single follow-up tuning settles an ambiguous real source
  atlas       synthetic error-rate atlas for six ALMA setups
  traps       three published octave-type misidentifications, rebuilt from their coverage
  tuning      tuning optimiser for ALMA Band 2 (now and with the Wideband Sensitivity Upgrade)
  ambiguity   the undertone family of candidate redshifts, the ambiguity function, and a sonification

Run any of them as  python -m ladder_redshifts.<module>  from the repository root.
"""
__version__ = '0.1.0'
