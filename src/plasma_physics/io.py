"""Load the processed TRICE-2 High data shipped in data/processed/.

Time is in seconds on the flight clock used by the original files.
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import scipy.io as sio

DATA_DIR = Path(__file__).resolve().parents[2] / "data" / "processed"


@dataclass(frozen=True)
class TimeSeries:
    """A uniformly sampled signal: timestamps `t` (s) and values `x`."""

    t: np.ndarray
    x: np.ndarray

    @property
    def fs(self) -> float:
        """Sampling rate (Hz), from the median time step."""
        return 1.0 / float(np.median(np.diff(self.t)))


def load_density(data_dir: Path = DATA_DIR) -> TimeSeries:
    """Electron density n_e (cm^-3), derived from the traced plasma frequency.

    Columns of chargeProfile: time (s), plasma frequency (kHz), n_e = (f_pe / 8.98)^2.
    """
    cp = sio.loadmat(data_dir / "charge_profile.mat")["chargeProfile"]
    return TimeSeries(t=cp[:, 0], x=cp[:, 2])


def load_efield(data_dir: Path = DATA_DIR) -> TimeSeries:
    """Electric field: one channel, apparently raw instrument counts (geometry undocumented)."""
    raw = sio.loadmat(data_dir / "EF_RawData.mat")["raw_data"]
    return TimeSeries(t=raw[:, 0], x=raw[:, 1])


def load_legacy_slopes(data_dir: Path = DATA_DIR) -> np.ndarray:
    """All 165 windows computed in MATLAB (slopes.txt).

    Columns: density slope, E-field slope, window start (s), window end (s).
    """
    return np.loadtxt(data_dir / "slopes.txt")


def load_selected_slopes(data_dir: Path = DATA_DIR) -> np.ndarray:
    """The 64 windows kept after manual inspection (my_slopes.mat), same columns."""
    return sio.loadmat(data_dir / "my_slopes.mat")["my_slopes"]
