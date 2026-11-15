import numpy as np

from plasma_physics.stats import block_bootstrap_ci


def test_block_bootstrap_ci_covers_the_true_mean_and_widens_with_autocorrelation():
    rng = np.random.default_rng(0)
    white = rng.normal(1.0, 1.0, 400)
    lo, hi = block_bootstrap_ci(white, block=1, n_boot=2000, rng=rng)
    assert lo < 1.0 < hi

    # AR(1) series: a block bootstrap must give a wider interval than an iid one
    ar = np.zeros(400)
    for i in range(1, 400):
        ar[i] = 0.8 * ar[i - 1] + rng.normal()
    iid = block_bootstrap_ci(ar, block=1, n_boot=2000, rng=rng)
    blocked = block_bootstrap_ci(ar, block=10, n_boot=2000, rng=rng)
    assert (blocked[1] - blocked[0]) > 1.5 * (iid[1] - iid[0])
