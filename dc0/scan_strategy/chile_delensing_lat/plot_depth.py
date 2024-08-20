import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np

# delens  priority_1000  priority_10000000  wide
m0 = hp.read_map("outputs/wide/mapmaker_cov.fits")
m1 = hp.read_map("outputs/priority_1000/mapmaker_cov.fits")
m2 = hp.read_map("outputs/priority_10000000/mapmaker_cov.fits")
m3 = hp.read_map("outputs/delens/mapmaker_cov.fits")

maps = np.vstack([m0, m1, m2, m3])

for m in maps:
    good = m != 0
    m[good] = 1 / np.sqrt(m[good])
    # m[good] /= np.amax(m[good])
    m[m == 0] = hp.UNSEEN

norm = np.amax(maps)
for m in maps:
    m[m != hp.UNSEEN] /= norm

nrow, ncol = 2, 2
fig = plt.figure(figsize=[6 * ncol, 4 * nrow])
for i, m in enumerate(maps):
    title = [
        "Wide",
        "Mostly Wide",
        "Mostly Delens",
        "Delens",
    ][i]
    hp.mollview(
        m,
        sub=[nrow, ncol, 1 + i],
        cmap="magma",
        title=title,
        # cbar=False,
        min=0,
        max=2,
    )

plt.savefig("relative_depths_delensing_chlat.png")

nrow, ncol = 1, 2
fig = plt.figure(figsize=[6 * ncol, 4 * nrow])
for i, m in enumerate([np.sqrt(2) * maps[0], 2 * maps[1]]):
    title = [
        "2 x Wide",
        "4 x Mostly Wide",
    ][i]
    hp.mollview(
        m,
        sub=[nrow, ncol, 1 + i],
        cmap="magma",
        title=title,
        # cbar=False,
        min=0,
        max=2,
    )

plt.savefig("relative_depths_hybrid_vs_wide.png")
