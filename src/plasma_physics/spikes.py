"""Automatic detection of the "sun spikes" in the E field.

The spikes are short rectangular pulses (~40 ms) every ~0.92 s, half the period of the large
~1.84 s oscillation (most likely the rocket spin), alternating in sign. In 2022 the clean
stretches between them were found by hand (11 starting points). Here every pulse is found from
its sharp edges.
"""

from dataclasses import dataclass

import numpy as np

from plasma_physics.io import TimeSeries


@dataclass(frozen=True)
class Spikes:
    start: np.ndarray  # time of the first sharp edge of each pulse (s)
    end: np.ndarray  # time of the last sharp edge (s)
    polarity: np.ndarray  # +1 pulse up, -1 pulse down, 0 undetermined


@dataclass(frozen=True)
class Segments:
    start: np.ndarray
    end: np.ndarray
    polarity: np.ndarray  # polarity of the spike right before the segment


def detect_spikes(efield: TimeSeries, k_mad: float = 12.0, merge_gap: float = 0.15) -> Spikes:
    """Find pulses as clusters of sample-to-sample jumps far above the typical jump.

    A jump is "sharp" when |dx| > median + k_mad * (robust std of |dx|). Edges closer than
    `merge_gap` seconds belong to the same pulse (its rising and falling edge).
    """
    jumps = np.abs(np.diff(efield.x))
    robust_std = 1.4826 * np.median(np.abs(jumps - np.median(jumps)))
    edges = efield.t[1:][jumps > np.median(jumps) + k_mad * robust_std]
    if edges.size == 0:
        empty = np.array([])
        return Spikes(empty, empty, empty)
    groups = np.split(edges, np.flatnonzero(np.diff(edges) > merge_gap) + 1)
    start = np.array([g[0] for g in groups])
    end = np.array([g[-1] for g in groups])
    return Spikes(start, end, _polarity(efield, start, end))


def _polarity(efield: TimeSeries, start: np.ndarray, end: np.ndarray) -> np.ndarray:
    """Sign of the pulse: median inside it minus the median of the 5-50 ms before it."""
    out = np.zeros(len(start))
    for k, (a, b) in enumerate(zip(start, end, strict=True)):
        inside = efield.x[(efield.t >= a) & (efield.t <= b)]
        before = efield.x[(efield.t >= a - 0.05) & (efield.t < a - 0.005)]
        if inside.size and before.size:
            out[k] = np.sign(np.median(inside) - np.median(before))
    return out


def clean_segments(
    spikes: Spikes, margin: float = 0.03, min_len: float = 0.75, max_len: float = 0.9
) -> Segments:
    """Stretches between consecutive spikes, trimmed by `margin` on both sides.

    Keeping only lengths in [min_len, max_len] drops gaps where a spike was missed
    (they are ~1.8 s long) and false detections (too short).
    """
    start = spikes.end[:-1] + margin
    end = spikes.start[1:] - margin
    length = end - start
    keep = (length >= min_len) & (length <= max_len)
    return Segments(start[keep], end[keep], spikes.polarity[:-1][keep])
