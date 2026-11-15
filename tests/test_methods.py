"""Synthetic checks: each method must recover what we put in."""

import numpy as np
import pytest

from plasma_physics.coherence import (
    coherence,
    coherency_test,
    permutation_test,
    phase_randomization_test,
    segment_coherency,
    segment_spectra,
)
from plasma_physics.io import TimeSeries
from plasma_physics.spectral import hann_periodogram, loglog_fit
from plasma_physics.spikes import Segments, clean_segments, detect_spikes


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


def unrelated_spectra_with_shared_bursts(rng: np.random.Generator):
    """Two signals with NO phase relation whose power rises together in a few segments."""
    n_seg, n_freq = 300, 20
    gain = np.where(rng.random(n_seg) < 0.05, 30.0, 1.0)[:, None]  # ~5% strong segments
    X = gain * (rng.normal(size=(n_seg, n_freq)) + 1j * rng.normal(size=(n_seg, n_freq)))
    Y = gain * (rng.normal(size=(n_seg, n_freq)) + 1j * rng.normal(size=(n_seg, n_freq)))
    return X, Y, np.ones(n_freq, dtype=bool)


def test_shuffle_null_is_fooled_by_shared_power_bursts_but_the_other_tests_are_not():
    """The reason the first analysis overstated significance, as a regression test."""
    shuffle_hits, phase_hits, coherency_hits = 0, 0, 0
    n_trials = 40
    for seed in range(n_trials):
        rng = np.random.default_rng(seed)
        X, Y, band = unrelated_spectra_with_shared_bursts(rng)
        shuffle_hits += permutation_test(X, Y, band, 200, rng)[1] < 0.05
        phase_hits += phase_randomization_test(X, Y, band, 200, rng)[1] < 0.05
        result = coherency_test(X, Y, band, 200, rng)
        coherency_hits += min(result.p_real, result.p_imag) < 0.025
    assert shuffle_hits > n_trials * 0.5  # far above the nominal 5% false-positive rate
    assert phase_hits <= n_trials * 0.2
    assert coherency_hits <= n_trials * 0.2


def test_segment_coherency_separates_in_phase_from_quadrature():
    rng = np.random.default_rng(5)
    X = rng.normal(size=(100, 20)) + 1j * rng.normal(size=(100, 20))
    band = np.ones(20, dtype=bool)
    assert np.allclose(segment_coherency(X, 2 * X, band), 1.0)
    assert np.allclose(segment_coherency(X, -1j * X, band), 1j)  # Y lags X by 90 degrees


def test_segment_spectra_recover_quadrature_through_a_flipping_sensor():
    """End to end: E is the time derivative of the density signal, seen by a sensor whose
    sign flips from one segment to the next. Only the sign correction recovers it."""
    rng = np.random.default_rng(6)
    fs_n, fs_e, total = 152.6, 2500.0, 60.0
    t_e = np.arange(0, total, 1 / fs_e)
    common = np.cumsum(rng.normal(size=t_e.size))  # red noise shared by both signals
    common -= np.linspace(common[0], common[-1], common.size)
    derivative = np.gradient(common, t_e)
    starts = np.arange(1.0, total - 2.0, 0.92)
    polarity = np.where(np.arange(len(starts)) % 2 == 0, 1.0, -1.0)
    sensor_sign = np.ones_like(t_e)
    for a, pol in zip(starts, polarity, strict=True):
        sensor_sign[(t_e >= a) & (t_e < a + 0.92)] = pol
    efield = TimeSeries(t_e, sensor_sign * derivative + rng.normal(0, 0.1, t_e.size))
    t_n = np.arange(0, total, 1 / fs_n)
    density = TimeSeries(t_n, 1000 + np.interp(t_n, t_e, common))
    segments = Segments(starts, starts + 0.8, polarity)

    def mean_coherency(sign_correct: bool) -> complex:
        spec = segment_spectra(density, efield, segments, sign_correct=sign_correct)
        band = (spec.f >= 10) & (spec.f <= 40)
        return segment_coherency(spec.X, spec.Y, band).mean()

    corrected = mean_coherency(True)
    assert corrected.imag < -0.5  # the derivative leads the signal by 90 degrees
    assert abs(corrected.real) < 0.2
    assert abs(mean_coherency(False)) < 0.2  # without the correction it cancels out
