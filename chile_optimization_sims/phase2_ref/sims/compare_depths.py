import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np


# f_quality from https://docs.google.com/spreadsheets/d/1vSUabA4JkrFM7xJDlhcrH5enWgxFbzuVKk-t0sswi-4/edit?usp=sharing

# Compare SAT maps

yield_ = 0.8

scale = {
    # f_quality * ntube * nyear
    "spsat" : {
        "f030" : 0.72 * 1 * 10,
        "f040" : 0.72 * 1 * 10,
        "f085" : 0.72 * 3 * 10,
        "f095" : 0.72 * 3 * 10,
        "f145" : 0.71 * 3 * 10,
        "f155" : 0.71 * 3 * 10,
        "f220" : 0.50 * 2 * 10,
        "f280" : 0.35 * 2 * 10,
    },
    "chsat" : {
        "f030" : 0.72 * 3 * 10,
        "f040" : 0.72 * 3 * 10,
        "f085" : 0.72 * 8 * 10,
        "f095" : 0.72 * 8 * 10,
        "f145" : 0.71 * 8 * 10,
        "f155" : 0.71 * 8 * 10,
        "f220" : 0.50 * 8 * 10,
        "f280" : 0.35 * 8 * 10,
    },
    # f_quality * ntele * nyear
    "splat" : {
        "f020" : 0.74 * 1 * 10,
        "f030" : 0.74 * 1 * 10,
        "f040" : 0.74 * 1 * 10,
        "f090" : 0.74 * 1 * 10,
        "f150" : 0.69 * 1 * 10,
        "f220" : 0.68 * 1 * 10,
        "f280" : 0.47 * 1 * 10,
    },
    "chlat" : {
        "f020" : 0.80 * 36,
        "f030" : 0.80 * 36,
        "f040" : 0.80 * 36,
        "f090" : 0.80 * 36,
        "f150" : 0.80 * 36,
        "f220" : 0.79 * 36,
        "f280" : 0.54 * 36,
    },
}


def invert_map(m):
    good = m != 0
    minv = np.zeros_like(m)
    minv[good] = 1 / m[good]
    return minv


def load_map(path):
    fname_season = f"{path}/season/mapmaker_cov.fits"
    fname_break = f"{path}/break/mapmaker_cov.fits"
    print(f"Loading {fname_season}")
    mseason = hp.read_map(fname_season)
    print(f"Loading {fname_break}")
    mbreak = hp.read_map(fname_break)
    return invert_map(invert_map(mseason) + invert_map(mbreak))


def plot_map(m, title, iplot, ref=None):
    # Transform pixel variance to depth
    mm = m.copy()
    nside = hp.get_nside(mm)
    area = hp.nside2pixarea(nside, degrees=True)
    mm = np.sqrt(mm * area) * 1e6 * 60  # uK.arcmin
    mm *= np.sqrt(2)  # I -> P
    good = m != 0
    mm[m == 0] = hp.UNSEEN
    mmsorted = np.sort(mm[good])
    limit = mmsorted[int(0.03 * mm.size)]
    best = np.logical_and(good, mm < limit)
    avg = np.mean(mm[best])
    title = title + f" {avg:.3f}"
    if ref is not None:
        title = title + f" (time : {(ref / avg)**2:.3f})"
    vmin = np.amin(mm[good])
    hp.mollview(
        mm,
        min=vmin,
        max=2 * vmin,
        cmap="magma",
        title=title,
        sub=[nrow, ncol, iplot],
        unit=r"$\mu$K.arcmin",
        format="%.3f",
    )
    return avg


bands = ["f030", "f040", "f085", "f095", "f145", "f155", "f220", "f280"]
nband = len(bands)

nrow, ncol = 2, nband
fig = plt.figure(figsize=[4 * ncol, 4 * nrow])
for iband, band in enumerate(bands):
    m0 = load_map(
        f"outputs/sun90/{band}"
    ) / scale["spsat"][band] / yield_
    m1 = load_map(
        f"../../phase2/sims/outputs/sun90/{band}"
    ) / scale["chsat"][band] / yield_
    avg = plot_map(m0, f"SPSAT {band}", iband + 1)
    plot_map(m1, f"CHSAT {band}", iband + 1 + ncol, avg)

fig.savefig("sat_depths.png")


bands = ["f020", "f030", "f040", "f090", "f150", "f220", "f280"]
nband = len(bands)

nrow, ncol = 2, nband
fig = plt.figure(figsize=[4 * ncol, 4 * nrow])
for iband, band in enumerate(bands):
    m0 = load_map(
        f"outputs/splat/{band}"
    ) / scale["splat"][band] / yield_
    m1 = load_map(
        f"../../phase2/sims/outputs/lat_delensing_tiled/{band}"
    ) / scale["chlat"][band] / yield_
    avg = plot_map(m0, f"SPLAT {band}", iband + 1)
    plot_map(m1, f"CHLAT {band}", iband + 1 + ncol, avg)

fig.savefig("lat_depths.png")
