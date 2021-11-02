import os
import sys

import healpy as hp
import numpy as np
import matplotlib.pyplot as plt

#rootdir = "/global/cscratch1/sd/zonca/cmbs4/map_based_simulations/202102_design_tool_run"
rootdir = "/global/cfs/cdirs/cmbs4/dm/dstool_202102/output/"
flavors = "noise", "foregrounds", "atmo", 

sites = [
    "pole", "pole", "pole", "chile",
]

telescopes = [
    "SAT", "SAT", "LAT", "LAT",
]

bands = [
    [None,     None,  "ULFPL1",    None],
    ["LFS1",   None,   "LFPL1",  "LFL1"],
    ["LFS2",   None,   "LFPL2",  "LFL2"],
    ["MFLS1", "MFHS1", "MFPL1",  "MFL1"],
    ["MFLS2", "MFHS2", "MFPL2",  "MFL2"],
    ["HFS1",   None,   "HFPL1",  "HFL1"],
    ["HFS2",   None,   "HFPL2",  "HFL2"],
]

freqs = [
    [None,  None,  20, None],
    [  27,  None,  26,   26],
    [  39,  None,  39,   39],
    [  85,    95,  92,   92],
    [ 145,   155, 149,  149],
    [ 225, None,  227,  227],
    [ 278, None,  286,  286],
]

nest = False
xsize = 1600

ncol = len(sites)
nrow = len(bands)

fig = plt.figure(figsize=[ncol * 4, nrow * 2])

for row in range(nrow):
    for col in range(ncol):
        #if col > 1:  # DEBUG
        #    continue
        site = sites[col]
        tele = telescopes[col]
        band = bands[row][col]
        freq = freqs[row][col]
        if band is None:
            continue
        if tele == "SAT":
            nside = 512
            vmax = 0.5
            if freq < 50 or freq > 200:
                vmax *= 4
        else:
            nside = 4096
            if site == "pole":
                vmax = 2
            else:
                vmax = 10
            if freq < 50 or freq > 200:
                vmax *= 5
        mtot = None
        for flavor in flavors:
            #if flavor in ["atmo"]:  # DEBUG
            #    continue
            fname = os.path.join(
                rootdir,
                flavor,
                f"{tele}-{band}_{site}",
                f"cmbs4_KCMB_{tele}-{band}_{site}_nside{nside}_1_of_1.fits",
            )
            print(fname)
            m = hp.read_map(fname, [1, 2], nest=nest)
            if mtot is None:
                mtot = m
            else:
                mtot += m
        if site == "pole":
            title1 = f"SP{tele}"
        else:
            title1 = f"CH{tele}"
        title2 = f"{freq}GHz"
        bad = np.abs(mtot[0]) > 1e10
        p = np.sqrt(np.sum(mtot ** 2, 0)) * 1e6
        p[bad] = hp.UNSEEN
        hp.mollview(
            p,
            sub=[nrow, ncol, row * ncol + col + 1],
            title=None,
            nest=nest,
            min=0,
            max=vmax,
            unit="$\mu$K$_\mathrm{CMB}$",
            xsize=xsize,
            cmap="inferno",
        )
        ax = plt.gca()
        ax.text(-0.1, 1.1, title1, horizontalalignment="left", verticalalignment="top", transform=ax.transAxes, size=18)
        ax.text(1.15, 1.1, title2, horizontalalignment="right", verticalalignment="top", transform=ax.transAxes, size=18)

fig.subplots_adjust(wspace=0.05, hspace=0.05, bottom=0.35, top=0.95)
fig.savefig("all_bands.pdf")
fig.savefig("all_bands.png")
