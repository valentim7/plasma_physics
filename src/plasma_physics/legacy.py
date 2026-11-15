"""Faithful port of matlab/step06_spectral_slopes.m (2022).

Kept deliberately close to the MATLAB code, quirks included, so that its output can be
checked number by number against data/processed/slopes.txt. Improved methods live in
the other modules.
"""

import numpy as np

from plasma_physics.io import TimeSeries

# Starting times picked by hand, each right after a sun spike (see docs/pipeline.md, step 5).
MANUAL_STARTS: tuple[float, ...] = (
    500.75,
    531.17,
    559.76,
    587.43,
    615.97,
    644.6,
    673.12,
    701.71,
    730.35,
    758.87,
    787.47,
)
DENSITY_FS = 1 / 0.0066  # MATLAB value; the true step is 0.006554 s (0.7% off, kept on purpose)
EFIELD_FS = 1 / 0.0004  # 2500 Hz
FIT_LOG10_RANGE = (1.0, 1.6)  # fit between 10 Hz and ~40 Hz


def matlab_hanning(n: int) -> np.ndarray:
    """MATLAB's hanning(n): symmetric Hann window without the zero end points."""
    k = np.arange(1, n + 1)
    return 0.5 * (1 - np.cos(2 * np.pi * k / (n + 1)))


def legacy_intervals(
    start: float, n_intervals: int = 30, step: float = 0.92, length: float = 0.75
) -> np.ndarray:
    """The 30 E-field sub-windows (start, end) used by the MATLAB script.

    Quirk reproduced: window 1 starts at `start`, but window i >= 2 starts at
    start + step * i (not step * (i - 1)), so there is a 1.84 s jump after window 1.
    """
    starts = np.array([start] + [step * i + start for i in range(2, n_intervals + 1)])
    return np.column_stack([starts, starts + length])


def _legacy_slope(log_f: np.ndarray, log_p: np.ndarray) -> float:
    """Least-squares slope from the first point with log10 f >= 1 to the first >= 1.6."""
    lo, hi = FIT_LOG10_RANGE
    i0 = int(np.argmax(log_f >= lo))
    i1 = int(np.argmax(log_f >= hi))
    return float(np.polyfit(log_f[i0 : i1 + 1], log_p[i0 : i1 + 1], 1)[0])


def legacy_slopes(start: float, density: TimeSeries, efield: TimeSeries) -> np.ndarray:
    """Return 15 rows of (density slope, E-field slope, window start, window end)."""
    intervals = legacy_intervals(start)

    # E field: same number of samples in every sub-window (the shortest one wins)
    idx = [np.flatnonzero((efield.t > a) & (efield.t < b)) for a, b in intervals]
    n_e = min(len(i) for i in idx)
    w = matlab_hanning(n_e)
    power_e = np.array([np.abs(np.fft.fft(w * efield.x[i[:n_e]])) ** 2 for i in idx]).T
    f_e = np.arange(n_e) * EFIELD_FS / n_e

    rows = []
    for r in range(len(intervals) // 2):
        # each output row pairs sub-windows 2r and 2r+1
        t0, t1 = intervals[2 * r, 0], intervals[2 * r + 1, 1]

        # density over the whole span of the pair, including the first sample >= t1
        i0 = int(np.argmax(density.t >= t0))
        i1 = int(np.argmax(density.t >= t1))
        d = density.x[i0 : i1 + 1]
        n = len(d)
        p = np.abs(np.fft.fft(matlab_hanning(n) * d)) ** 2
        f = np.arange(n) * DENSITY_FS / n
        top = n // 2
        slope_n = _legacy_slope(np.log10(f[1:top]), np.log10(p[1:top]))

        # E field: average the two sub-window spectra, then fit
        p_e = power_e[:, 2 * r : 2 * r + 2].mean(axis=1)
        top_e = n_e // 2
        slope_e = _legacy_slope(np.log10(f_e[1:top_e]), np.log10(p_e[1:top_e]))

        rows.append((slope_n, slope_e, t0, t1))
    return np.array(rows)


def legacy_all(density: TimeSeries, efield: TimeSeries) -> np.ndarray:
    """All 165 rows for the 11 manual starting points."""
    return np.vstack([legacy_slopes(s, density, efield) for s in MANUAL_STARTS])
