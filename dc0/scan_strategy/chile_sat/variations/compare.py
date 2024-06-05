import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np


flavors = "baseline", "no_break", "no_lunar_avoidance", "south_only", "el_min_40", "sidereal", "all"

depths = {}
for flavor in flavors:
    try:
        cov = hp.read_map(f"outputs/{flavor}/mapmaker_cov.fits")
    except:
        cov = hp.read_map(f"outputs_old/{flavor}/mapmaker_cov.fits")
    depth = np.sqrt(cov)
    depths[flavor] = depth

hits = {}
for flavor in flavors:
    try:
        hits[flavor] = hp.read_map(f"outputs/{flavor}/mapmaker_hits.fits")
    except:
        hits[flavor] = hp.read_map(f"outputs_old/{flavor}/mapmaker_hits.fits")

nrow, ncol = 2, 4
fig = plt.figure(figsize=[4 * ncol, 4 * nrow])
for i, flavor in enumerate(flavors):
    baseline = depths["baseline"]
    if flavor == "baseline":
        depth = baseline.copy()
        good = depth != 0
        depth[good] /= np.amin(depth[good])
        depth[np.logical_not(good)] = hp.UNSEEN
        hp.mollview(depth, min=1, max=4, title=f"{flavor} / min(baseline)", sub=[nrow, ncol, 1 + i], cmap="magma")
    else:
        hit_ratio = np.sum(hits[flavor]) / np.sum(hits["baseline"])
        ref = depths[flavor]
        ratio = np.zeros_like(baseline)
        good = np.logical_and(baseline != 0, ref != 0)
        ratio[good] = ref[good] / baseline[good]
        ratio[np.logical_not(good)] = hp.UNSEEN
        hp.mollview(
            ratio,
            min=0.7,
            max=1.3,
            title=f"{flavor} / baseline, {hit_ratio:.3f}",
            sub=[nrow, ncol, 1 + i],
            cmap="bwr",
        )

fig.savefig("relative_depth.png")
plt.show()
