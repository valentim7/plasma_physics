"""Power spectra and power-law (log-log) slope fits."""

from dataclasses import dataclass

import numpy as np
from scipy.signal import detrend as _detrend

FIT_BAND_HZ = (10.0, 10**1.6)  # same band as the 2022 analysis: 10 Hz to ~40 Hz


@dataclass(frozen=True)
class SlopeFit:
    slope: float
    stderr: float  # OLS standard error; assumes independent log-power residuals
    n_points: int


def hann_periodogram(
    x: np.ndarray, fs: float, detrend: bool = True
) -> tuple[np.ndarray, np.ndarray]:
    """One-sided Hann-windowed periodogram (unnormalized power, enough for slopes).

    `detrend` removes a straight line first, so the mean and slow drifts (for example the
    rocket spin in the E field) do not leak into the fit band.
    """
    x = _detrend(x, type="linear") if detrend else np.asarray(x, dtype=float)
    w = np.hanning(len(x))
    p = np.abs(np.fft.rfft(w * x)) ** 2
    f = np.fft.rfftfreq(len(x), 1 / fs)
    return f, p


def loglog_fit(f: np.ndarray, p: np.ndarray, band: tuple[float, float] = FIT_BAND_HZ) -> SlopeFit:
    """Fit log10 p = a + slope * log10 f inside `band` (inclusive) by least squares."""
    sel = (f >= band[0]) & (f <= band[1]) & (p > 0)
    x, y = np.log10(f[sel]), np.log10(p[sel])
    (slope, intercept), cov = np.polyfit(x, y, 1, cov=True)
    return SlopeFit(slope=float(slope), stderr=float(np.sqrt(cov[0, 0])), n_points=int(sel.sum()))
