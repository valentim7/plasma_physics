"""Cross-spectral coherence between density and E field over clean segments.

If density and potential follow the Boltzmann relation, density fluctuations and E-field
fluctuations are linearly related, so their coherence should be clearly above zero. A slope
comparison alone cannot show that; coherence can.
"""

from dataclasses import dataclass

import numpy as np
from scipy.signal import butter, detrend, sosfiltfilt

from plasma_physics.io import TimeSeries
from plasma_physics.spikes import Segments


@dataclass(frozen=True)
class SegmentSpectra:
    f: np.ndarray  # frequencies (Hz)
    X: np.ndarray  # density FFTs, one row per segment
    Y: np.ndarray  # E-field FFTs (sign-corrected), one row per segment
    t: np.ndarray  # segment start times (s)


def segment_spectra(
    density: TimeSeries,
    efield: TimeSeries,
    segments: Segments,
    duration: float = 0.75,
    lowpass_hz: float = 60.0,
    sign_correct: bool = True,
) -> SegmentSpectra:
    """FFT both signals on the density time grid, segment by segment.

    The E field (2500 Hz) is low-pass filtered (zero phase) below the density Nyquist
    (~76 Hz) and interpolated onto the density timestamps. With `sign_correct`, each E
    segment is multiplied by the polarity of the preceding spike. The spikes alternate in sign,
    which suggests the sensor flips orientation every half spin; without the flip the
    cross-spectrum cancels. The instrument geometry is inferred from the signal.
    """
    sos = butter(6, lowpass_hz, fs=efield.fs, output="sos")
    e_smooth = sosfiltfilt(sos, efield.x - np.median(efield.x))
    e_on_density_grid = np.interp(density.t, efield.t, e_smooth)

    n = int(duration * density.fs)
    w = np.hanning(n)
    X, Y, t = [], [], []
    for a, pol in zip(segments.start, segments.polarity, strict=True):
        i0 = int(np.searchsorted(density.t, a))
        if a < density.t[0] or i0 + n > len(density.t) or a > efield.t[-1] - duration:
            continue
        x = density.x[i0 : i0 + n]
        x = detrend(x / x.mean())  # relative density fluctuation, delta n / n
        y = detrend(e_on_density_grid[i0 : i0 + n])
        if sign_correct:
            y = y * (pol if pol != 0 else 1.0)
        X.append(np.fft.rfft(w * x))
        Y.append(np.fft.rfft(w * y))
        t.append(a)
    f = np.fft.rfftfreq(n, 1 / density.fs)
    return SegmentSpectra(f, np.array(X), np.array(Y), np.array(t))


def coherence(X: np.ndarray, Y: np.ndarray) -> np.ndarray:
    """Magnitude-squared coherence per frequency, averaging over segments (rows)."""
    sxy = (X * np.conj(Y)).sum(axis=0)
    sxx = (np.abs(X) ** 2).sum(axis=0)
    syy = (np.abs(Y) ** 2).sum(axis=0)
    return np.abs(sxy) ** 2 / (sxx * syy)


def band_coherence(X: np.ndarray, Y: np.ndarray, band_mask: np.ndarray) -> float:
    return float(coherence(X, Y)[band_mask].mean())


def permutation_test(
    X: np.ndarray, Y: np.ndarray, band_mask: np.ndarray, n_perm: int, rng: np.random.Generator
) -> tuple[float, float, np.ndarray]:
    """Observed band coherence and its permutation p-value.

    The null pairs each density segment with a random E segment from the same set, which
    keeps both spectra intact but breaks any real time alignment.
    """
    observed = band_coherence(X, Y, band_mask)
    null = np.array(
        [band_coherence(X, Y[rng.permutation(len(Y))], band_mask) for _ in range(n_perm)]
    )
    p_value = (1 + np.sum(null >= observed)) / (1 + n_perm)
    return observed, float(p_value), null
