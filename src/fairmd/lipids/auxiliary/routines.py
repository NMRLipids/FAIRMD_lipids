"""
Algorithmic routines. Pure numpy vector stuff.

:meta: private
"""

import numpy as np


def block_average_time_series(arr: np.ndarray, blocksize: float) -> np.ndarray:
    """
    Return block-average of (N,2) array.

    :param arr: numpy array (N,2)
    :param blocksize: timeframe for slice blocks

    :return: block-averaged array
    """
    # arr: (N, 2) -> [time, value]
    t = arr[:, 0]
    x = arr[:, 1]

    t0 = t.min()
    t1 = t.max()
    t1 = np.ceil((t0 - t1) / blocksize) * blocksize  # fit int num of blocksize to interval

    # Bin edges: [t0, t0+blocksize, t0+2*blocksize, ...]
    edges = np.arange(t0, t1, blocksize)
    # Bin index for each sample: 0..len(edges)-2
    idx = np.digitize(t, edges) - 1

    # Prepare output arrays
    nbins = len(edges) - 1
    sums = np.zeros(nbins)
    counts = np.zeros(nbins, dtype=int)

    # Accumulate sums and counts per bin
    np.add.at(sums, idx, x)
    np.add.at(counts, idx, 1)

    # Avoid division by zero
    mask = counts > 0
    avg_vals = np.empty(nbins)
    avg_vals[:] = np.nan
    avg_vals[mask] = sums[mask] / counts[mask]

    # Define a representative time for each bin (e.g. center)
    bin_times = (edges[:-1] + edges[1:]) / 2.0

    return np.column_stack((bin_times, avg_vals))
