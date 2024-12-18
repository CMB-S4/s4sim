import glob
import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np


def invert_map(m):
    """Invert non-zero pixels"""
    result = np.zeros_like(m)
    good = m != 0
    result[good] = 1 / m[good]
    return result

lmax = 1024

bands = {
    "f026" : "f030",
    "f039" : "f040",
    "f085" : "f085",
    "f095" : "f095",
    "f145" : "f145",
    "f155" : "f155",
    "f227" : "f220",
    "f286" : "f280",
}

nrow, ncol = 2, 4
fig = plt.figure(figsize=[ncol * 4, nrow * 3])
iplot = 0
for alt_band, band in bands.items():
    iplot += 1
    print(f"{alt_band} {band}")
    fname_rhit = f"rhits/rhits_sat_{band}.fits"
    rhit = hp.read_map(fname_rhit)
    # Use best 3% for power spectrum
    limit = np.sort(rhit)[int(rhit.size * .97)]
    rhit[rhit < limit] = 0
    fsky = np.sum(rhit**2) / rhit.size

    levels0 = {}
    levels1 = {}
    cls0 = {}
    cls1 = {}

    nyear = 10
    fname1 = f"noise_{nyear:02}_years/phase2_noise_{alt_band}_SAT_mc_0000.fits"
    fname0 = "../phase1/" + fname1.replace("phase2", "phase1")
    alt_band = os.path.basename(fname1).split("_")[2]
    noise0 = hp.read_map(fname0, None)
    noise1 = hp.read_map(fname1, None)
    cl0 = hp.anafast(rhit * noise0, lmax=lmax, iter=0) / fsky
    cl1 = hp.anafast(rhit * noise1, lmax=lmax, iter=0) / fsky
    level0 = np.mean(cl0[2, 500:1000])
    level1 = np.mean(cl1[2, 500:1000])
    ratio = level1 / level0

    depth0 = np.sqrt(level0) * 1e6 * 180 / np.pi * 60
    depth1 = np.sqrt(level1) * 1e6 * 180 / np.pi * 60

    ax = fig.add_subplot(nrow, ncol, iplot)
    ax.set_title(f"{band} : " + r"C$_\ell^\mathrm{BB}$" + f" ratio = {ratio:.3f}")
    ax.loglog(cl0[2], label=f"Phase 1 : depth = {depth0:.3f}")
    ax.loglog(cl1[2], label=f"Phase 2 : depth = {depth1:.3f}")
    ax.legend(loc="best")
fig.tight_layout()
fig.savefig("sat_noise_comparison.png")
