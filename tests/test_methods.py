"""Synthetic checks: each method must recover what we put in."""

import numpy as np
import pytest

from plasma_physics.coherence import coherence, permutation_test
from plasma_physics.io import TimeSeries
from plasma_physics.spectral import hann_periodogram, loglog_fit
from plasma_physics.spikes import clean_segments, detect_spikes


def power_law_noise(n: int, fs: float, slope: float, rng: np.random.Generator) -> np.ndarray:
    """Gaussian noise whose power spectrum goes as f**slope."""
    f = np.fft.rfftfreq(n, 1 / fs)
    amp = np.zeros_like(f)
    amp[1:] = f[1:] ** (slope / 2)
    phases = np.exp(2j * np.pi * rng.random(len(f)))
    return np.fft.irfft(amp * phases * rng.normal(size=len(f)), n)


@pytest.mark.parametrize("true_slope", [-1.0, -2.0, -3.0])
def test_loglog_fit_recovers_known_power_law(true_slope):
    rng = np.random.default_rng(0)
    fs, n = 152.6, 4096
    spectra = [hann_periodogram(power_law_noise(n, fs, true_slope, rng), fs)[1] for _ in range(20)]
    f = np.fft.rfftfreq(n, 1 / fs)
    fit = loglog_fit(f, np.mean(spectra, axis=0))
    assert fit.slope == pytest.approx(true_slope, abs=0.1)
    assert fit.stderr > 0


def synthetic_efield(rng: np.random.Generator, n_pulses: int = 40) -> TimeSeries:
    """Spin sinusoid + noise + rectangular pulses every 0.92 s, alternating in sign."""
    fs = 2500.0
    t = np.arange(0, 0.92 * (n_pulses + 1), 1 / fs)
    x = 33000 + 1500 * np.sin(2 * np.pi * t / 1.84) + rng.normal(0, 5, t.size)
    for k in range(n_pulses):
        on = (t >= 0.5 + 0.92 * k) & (t < 0.54 + 0.92 * k)
        x[on] += 1500 * (1 if k % 2 == 0 else -1)
    return TimeSeries(t, x)


def test_detect_spikes_finds_every_pulse_and_its_sign():
    ef = synthetic_efield(np.random.default_rng(1))
    spikes = detect_spikes(ef)
    assert len(spikes.start) == 40
    assert np.allclose(np.diff(spikes.start), 0.92, atol=0.002)
    assert np.all(spikes.polarity[::2] == 1) and np.all(spikes.polarity[1::2] == -1)


def test_clean_segments_drop_a_gap_with_a_missed_spike():
    spikes = detect_spikes(synthetic_efield(np.random.default_rng(2)))
    full = clean_segments(spikes)
    assert len(full.start) == 39
    # removing one detection merges two gaps into a ~1.8 s gap, which must be rejected
    keep = np.ones(len(spikes.start), dtype=bool)
    keep[10] = False
    pruned = type(spikes)(spikes.start[keep], spikes.end[keep], spikes.polarity[keep])
    assert len(clean_segments(pruned).start) == 37


def test_coherence_is_one_for_linearly_related_and_low_for_independent_signals():
    rng = np.random.default_rng(3)
    X = rng.normal(size=(200, 30)) + 1j * rng.normal(size=(200, 30))
    assert np.allclose(coherence(X, -3 * X), 1.0)
    Y = rng.normal(size=(200, 30)) + 1j * rng.normal(size=(200, 30))
    assert coherence(X, Y).mean() < 0.02


def test_permutation_test_separates_signal_from_noise():
    rng = np.random.default_rng(4)
    X = rng.normal(size=(150, 20)) + 1j * rng.normal(size=(150, 20))
    noise = rng.normal(size=(150, 20)) + 1j * rng.normal(size=(150, 20))
    band = np.ones(20, dtype=bool)
    _, p_signal, _ = permutation_test(X, 0.4 * X + noise, band, 200, rng)
    _, p_noise, _ = permutation_test(X, noise, band, 200, rng)
    assert p_signal < 0.01
    assert p_noise > 0.05
