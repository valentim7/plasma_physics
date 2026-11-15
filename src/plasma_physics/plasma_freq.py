"""Next step (not implemented): trace the plasma frequency from raw HF spectra.

The HF spectrograms ("graydata", ~50M points for the flight) are not in this repository.
In 2021-2022 the cutoff was found with a hand-tuned threshold (the 0.997 factor in
matlab/step03_auto_plasma_frequency.m) and then corrected point by point with ginput
(step04). That manual step limits both the accuracy and the noise floor of the density.

What the existing hand-corrected trace tells us (data/processed/charge_profile.mat)
-------------------------------------------------------------------------------
The frequency bin is 0.3 kHz, but consecutive samples of the trace differ by 2 to 4 kHz
(median) and by up to ~13 kHz in 90% of cases; the tail reaches hundreds of kHz. Before
500 s the median step is ~1.4% of the plasma frequency, about 2.7% in density, which is
consistent with the flat density spectrum seen there at 10-40 Hz. The trace jitters by
many bins, so there is room to improve, but real fluctuations and tracing noise overlap in
the 15-50 kHz range.

Plan
----
0. Baseline first: a robust outlier filter (Hampel) on the 2022 semi-automatic trace. Any
   more complex method has to beat it.
1. For each spectrum, score candidate cutoff frequencies by the drop in power across them
   (a change-point score on the smoothed log spectrum).
2. Track the cutoff through time with a Viterbi path, constrained below the upper-hybrid
   cutoff. The jump penalty must be flat for small moves and grow only for large ones:
   the science signal is the small, fast fluctuation of the cutoff, and a penalty on every
   move would smooth it away.
3. Allow "no detection": where no candidate stands out (for example, no waves near the
   cutoff), return NaN instead of inventing a point.
4. Refine each point with a local step fit (sub-bin) and attach an uncertainty from the
   forward-backward posterior.
5. In parallel, try a differential estimate: the shift that best aligns consecutive spectra.
   For fluctuations it may be less noisy than locating the edge in each spectrum.

Validation
----------
- Synthetic spectrograms with a known cutoff and injected 10-40 Hz fluctuations: the
  recovered fluctuation spectrum must match the injected one (checks that step 2 does not
  smooth the signal).
- Density noise floor at 10-40 Hz, before and after, on real data. This is the main
  criterion.
- Agreement with the hand-corrected trace on a held-out part of the flight, only to catch
  gross errors: that trace carries its own biases and cannot prove the new one is better.
- The density/E-field coupling is computed only at the end and is never used for tuning.

Caveat: this plan was written without access to the HF spectra. The cutoff may not always
look like a clean edge, and the plan should be revised once the spectra are inspected.
"""

import numpy as np


def track_plasma_frequency(
    freq_khz: np.ndarray,
    log_power: np.ndarray,
    free_move_khz: float,
    jump_cost: float,
) -> np.ndarray:
    """Return one plasma-frequency estimate (kHz) per spectrum (column of `log_power`).

    Moves smaller than `free_move_khz` between consecutive spectra cost nothing; larger ones
    cost `jump_cost` per kHz beyond it. Spectra without a clear cutoff give NaN.
    """
    raise NotImplementedError("Needs the HF graydata; see the module docstring.")
