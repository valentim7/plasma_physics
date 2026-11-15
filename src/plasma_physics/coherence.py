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
    """Observed band coherence and its p-value under shuffled segment pairings.

    Caveat: this null is too lenient when the power of both signals rises and falls together
    (here, inside the aurora). Coherence weights segments by power, so "strong with strong"
    pairs inflate the observed value even with no phase relation, and shuffling removes that.
    Kept for comparison; use `phase_randomization_test` or `quadrature_test` instead.
    """
    observed = band_coherence(X, Y, band_mask)
    null = np.array(
        [band_coherence(X, Y[rng.permutation(len(Y))], band_mask) for _ in range(n_perm)]
    )
    p_value = (1 + np.sum(null >= observed)) / (1 + n_perm)
    return observed, float(p_value), null


def phase_randomization_test(
    X: np.ndarray, Y: np.ndarray, band_mask: np.ndarray, n_perm: int, rng: np.random.Generator
) -> tuple[float, float, np.ndarray]:
    """Band coherence tested against a null that keeps the pairing and the power of every
    segment, and only rotates each E segment by one random phase.

    Power co-variation is therefore the same in the null and in the data; only a phase
    relation that repeats from segment to segment can make the observed value stand out.
    """
    observed = band_coherence(X, Y, band_mask)
    null = np.empty(n_perm)
    for k in range(n_perm):
        rotation = np.exp(2j * np.pi * rng.random((len(Y), 1)))
        null[k] = band_coherence(X, Y * rotation, band_mask)
    p_value = (1 + np.sum(null >= observed)) / (1 + n_perm)
    return observed, float(p_value), null


def effective_segments(X: np.ndarray, Y: np.ndarray, band_mask: np.ndarray) -> float:
    """How many segments effectively carry the power-weighted coherence estimate.

    Equals the number of segments when all have the same power, and approaches 1 when a
    single segment dominates.
    """
    w = (np.abs(X[:, band_mask]) ** 2).mean(axis=1) * (np.abs(Y[:, band_mask]) ** 2).mean(axis=1)
    return float(w.sum() ** 2 / (w**2).sum())


def segment_coherency(X: np.ndarray, Y: np.ndarray, band_mask: np.ndarray) -> np.ndarray:
    """One complex number per segment: its band cross-spectrum divided by its own power.

    Every segment weighs the same, whatever its power. The real part measures in-phase
    coupling, the imaginary part coupling with a 90 degree phase shift (quadrature). With
    this sign convention (X times conj(Y)), a positive imaginary part means X leads Y.
    """
    cross = (X * np.conj(Y))[:, band_mask].sum(axis=1)
    px = (np.abs(X[:, band_mask]) ** 2).sum(axis=1)
    py = (np.abs(Y[:, band_mask]) ** 2).sum(axis=1)
    return cross / np.sqrt(px * py)


@dataclass(frozen=True)
class CoherencyTest:
    mean: complex  # mean per-segment coherency
    null_sd_real: float  # spread of the mean under shuffled pairings
    null_sd_imag: float
    p_real: float  # two-sided
    p_imag: float


def coherency_test(
    X: np.ndarray, Y: np.ndarray, band_mask: np.ndarray, n_perm: int, rng: np.random.Generator
) -> CoherencyTest:
    """Mean per-segment coherency, tested by shuffling the segment pairings.

    Because each pair is normalized by its own power, the power co-variation that breaks
    `permutation_test` does not affect this one. Shuffling also keeps any waveform that is
    locked to the segment start in both signals, so such an artifact would not look
    significant here.
    """
    observed = segment_coherency(X, Y, band_mask).mean()
    null = np.array(
        [segment_coherency(X, Y[rng.permutation(len(Y))], band_mask).mean() for _ in range(n_perm)]
    )
    p_real = (1 + np.sum(np.abs(null.real) >= abs(observed.real))) / (1 + n_perm)
    p_imag = (1 + np.sum(np.abs(null.imag) >= abs(observed.imag))) / (1 + n_perm)
    return CoherencyTest(
        mean=complex(observed),
        null_sd_real=float(null.real.std()),
        null_sd_imag=float(null.imag.std()),
        p_real=float(p_real),
        p_imag=float(p_imag),
    )
