"""Next step (not implemented): trace the plasma frequency from raw HF spectra.

The HF spectrograms ("graydata", ~50M points for the flight) are not in this repository.
In 2021-2022 the cutoff was found with a hand-tuned threshold (the 0.997 factor in
matlab/step03_auto_plasma_frequency.m) and then corrected point by point with ginput
(step04). That manual step limits both the accuracy and the noise floor of the density.

Planned approach once the HF data is available:
1. For each spectrum, score candidate cutoff frequencies by the drop in power across them
   (a change-point score on the smoothed log spectrum).
2. Track the cutoff through time with a Viterbi path that penalizes large jumps between
   consecutive spectra, so isolated spikes and dips cannot pull the trace off.
3. Validate against the hand-corrected trace in data/processed/charge_profile.mat
   (median absolute error in kHz, and the density noise floor at 10-40 Hz).
"""

import numpy as np


def track_plasma_frequency(
    freq_khz: np.ndarray, log_power: np.ndarray, max_jump_khz: float = 50.0
) -> np.ndarray:
    """Return one plasma-frequency estimate (kHz) per spectrum (column of `log_power`)."""
    raise NotImplementedError("Needs the HF graydata; see the module docstring.")
