import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np


cov1 = hp.read_map("../../phase1/sims/outputs/lat_wide/f090/mapmaker_cov.fits")
cov2 = hp.read_map("outputs/lat_wide/f090/mapmaker_cov.fits")

thinfp = 16
nyear = 16
area = hp.nside2pixarea(hp.get_nside(cov1), degrees=True)

depth1 = np.sqrt(cov1 * area  / thinfp / nyear) * 1e6 * 60
depth2 = np.sqrt(cov2 * area  / thinfp / nyear) * 1e6 * 60


nrow, ncol = 1, 2
fig = plt.figure(figsize=[ncol * 4, nrow * 4])

vmin = np.amin(depth1[depth1 != 0])
fsky = np.sum(depth1 != 0) / depth1.size
hp.mollview(
    depth1,
    min=vmin,
    max=2*vmin,
    sub=[nrow, ncol, 1],
    title=f"Phase 1, fsky={fsky:.3f}",
    format="%.3f",
    unit=f"$\mu$K.arcmin",
)

vmin = np.amin(depth2[depth2 != 0])
fsky = np.sum(depth2 != 0) / depth1.size
hp.mollview(
    depth2,
    min=vmin,
    max=2*vmin,
    sub=[nrow, ncol, 2],
    title=f"Phase 2, fsky={fsky:.3f}",
    format="%.3f",
    unit=f"$\mu$K.arcmin",
)

fig.savefig("plots/phase1_vs_phase2_depth.png")
