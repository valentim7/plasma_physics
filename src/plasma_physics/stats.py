"""Small statistics helpers."""

import numpy as np


def block_bootstrap_ci(
    x: np.ndarray,
    block: int,
    n_boot: int,
    rng: np.random.Generator,
    level: float = 0.95,
) -> tuple[float, float]:
    """Confidence interval for the mean with a moving-block bootstrap.

    Neighbouring windows are correlated in time, so resampling single points would make the
    interval too narrow; resampling runs of `block` consecutive points keeps that correlation.
    """
    x = np.asarray(x, dtype=float)
    n = len(x)
    n_blocks = int(np.ceil(n / block))
    starts = rng.integers(0, n - block + 1, size=(n_boot, n_blocks))
    idx = (starts[:, :, None] + np.arange(block)).reshape(n_boot, -1)[:, :n]
    means = x[idx].mean(axis=1)
    alpha = (1 - level) / 2
    lo, hi = np.quantile(means, [alpha, 1 - alpha])
    return float(lo), float(hi)
