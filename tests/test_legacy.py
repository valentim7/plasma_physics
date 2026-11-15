"""The Python port must reproduce the MATLAB results stored in slopes.txt."""

import numpy as np
import pytest
from scipy.signal.windows import hann

from plasma_physics.io import load_density, load_efield, load_legacy_slopes
from plasma_physics.legacy import legacy_all, legacy_intervals, matlab_hanning


def test_matlab_hanning_matches_symmetric_hann_without_zeros():
    assert np.allclose(matlab_hanning(9), hann(11)[1:-1])


def test_legacy_intervals_reproduce_the_jump_after_the_first_window():
    iv = legacy_intervals(500.75)
    assert iv.shape == (30, 2)
    assert iv[0, 0] == pytest.approx(500.75)
    assert iv[1, 0] == pytest.approx(500.75 + 1.84)
    assert np.allclose(iv[:, 1] - iv[:, 0], 0.75)


@pytest.fixture(scope="module")
def reproduced():
    return legacy_all(load_density(), load_efield())


def test_reproduces_every_matlab_window(reproduced):
    ref = load_legacy_slopes()

    def key(rows):
        return [(round(a, 2), round(b, 2)) for a, b in rows[:, 2:4]]

    mine = dict(zip(key(reproduced), reproduced[:, :2], strict=True))
    assert set(key(ref)) <= set(mine)
    diffs = np.array([mine[k] - r for k, r in zip(key(ref), ref[:, :2], strict=True)])
    # slopes.txt was written with %f (6 decimals), so agreement is limited by rounding
    assert np.abs(diffs).max() < 1e-5
