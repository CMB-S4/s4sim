import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np


m0 = hp.read_map("../chile_sat/variations/outputs/baseline/mapmaker_cov.fits")
m1 = hp.read_map("outputs/delens/mapmaker_cov.fits")
m2 = hp.read_map("core_only/outputs/mapmaker_cov.fits")
m3 = hp.read_map("full_footprint/outputs/mapmaker_cov.fits")

maps = np.vstack([m0, m1, m2, m3])

for m in maps:
    good = m != 0
    m[good] = 1 / np.sqrt(m[good])
    # m[good] /= np.amax(m[good])
    m[m == 0] = hp.UNSEEN

norm0 = np.amax(maps[0])
maps[0][maps[0] != hp.UNSEEN] /= norm0

norm = np.amax(maps[1:])
for m in maps[1:]:
    m[m != hp.UNSEEN] /= norm

nrow, ncol = 2, 2
fig = plt.figure(figsize=[6 * ncol, 4 * nrow])
for i, m in enumerate(maps):
    title = [
        "CHSAT",
        "Delens baseline",
        "Core only",
        "Full footprint",
    ][i]
    hp.mollview(
        m,
        sub=[nrow, ncol, 1 + i],
        cmap="magma",
        title=title,
        # cbar=False,
        min=0,
        max=1,
    )

plt.savefig("relative_depths_core_vs_full_chlat.png")
