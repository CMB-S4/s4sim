import numpy as np
import os
import sys
import matplotlib.pyplot as plt

import healpy as hp

# from toast_planck.utilities import plug_holes

plt.style.use("classic")

dec_min = -70
dec_max = 70
dec_step = 20
dec_overlap = 5
ra_step = 10
ra_overlap = 5

nbin = 10
# fwhm, tlim, plim = 10, 300 * 1e-6, 1 * 1e-7
fwhm, tlim, plim = 1, 500 * 1e-6, 2 * 1e-7
lmax = 512

wpatch = 10
hpatch = 20

south_tier_ranges = [(60, 30), (90, 10), (100, -65)]  # Tier 1, 2, 3
north_tier_ranges = [(165, 135), (175, 115), (230, 110)]  # Tier 1, 2, 3


def get_tier(left, right):
    tier = 4
    for ranges in [south_tier_ranges, north_tier_ranges]:
        for i, lims in enumerate(ranges):
            if left <= lims[0] and right >= lims[1]:
                return i + 1
    return tier


south_left = 100
south_top = [
    # Tier 3
    -15,  # 100
    -15,
    # Tier 2
    -20,  # 90
    -20,
    -25,  # 80
    # Tier 1
    -25,
    -30,  # 70
    -30,
    -30,  # 60
    # Tier 2
    -35,
    -35,  # 50
    -35,
    -30,  # 40
    -35,
    -35,  # 30
    -35,
    -35,  # 20
    # Tier 3
    -30,
    -30,  # 10
    -25,
    -25,  # 0,
    -25,
    -20,  # -10
    -20,
    -15,  # -20
    -15,
    -15,  # -30
    -10,
    -5,  # -40
    0,
    5,  # -50
    5,
]

north_left = 230
north_top = [
    # Tier 3
    -20,  # 230
    -25,
    -25,  # 220
    -25,
    -25,  # 210
    -25,
    -25,  # 200
    -20,
    -10,  # 190
    0,
    10,  # 180
    10,
    # Tier 2
    10,  # 170
    # Tier 1
    15,
    10,  # 160
    5,
    10,  # 150
    10,
    # Tier 2
    15,  # 140
    15,
    20,  # 130
    20,
    # Tier 3
    20,  # 120
]

nside = 512
npix = 12 * nside ** 2

hp.mollview(np.zeros(npix), cbar=False, title="Patch positions")
hp.graticule(22.5)

def plot_patch(patch, tier):
    left, top, right, bottom = patch
    color, lw = [("red", 3), ("black", 2), ("grey", 2)][tier - 1]
    lon = [left, left, right, right, left]
    lat = [bottom, top, top, bottom, bottom]
    n = len(lon)
    x = np.arange(n) / (n - 1)
    xfull = np.linspace(0, 1, 100 * n)
    lonfull = np.interp(xfull, x, lon)
    latfull = np.interp(xfull, x, lat)
    hp.projplot(
        lonfull,
        latfull,
        "-",
        threshold=1,
        lonlat=True,
        color=color,
        lw=lw,
        alpha=0.8,
        coord="C",
        zorder=(4 - tier) * 100,
    )


def add_patch(fout, patch, tier, mult=1):
    ra_stop, dec_start, ra_start, dec_stop = patch
    if ra_start < 0:
        ra_start += 360
    if ra_stop < 0:
        ra_stop += 360
    name = "Tier{}DEC{:+04}..{:+04}_RA{:+04}..{:+04}".format(
        tier, dec_start, dec_stop, ra_start, ra_stop
    )
    priority = [1, 1000, 1000000][tier - 1] * mult
    fout.write("--patch\n")
    fout.write(
        "{},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f}\n".format(
            name, priority, ra_start, dec_stop, ra_stop, dec_start
        )
    )


with open("patches.txt", "w") as fout:
    left = south_left
    for top in south_top:
        tier = get_tier(left, left - wpatch)
        patch = [left, top, left - wpatch, top - hpatch]
        plot_patch(patch, tier)
        add_patch(fout, patch, tier, mult=1)
        left -= 5

    left = north_left
    for top in north_top:
        tier = get_tier(left, left - wpatch)
        patch = [left, top, left - wpatch, top - hpatch]
        plot_patch(patch, tier)
        add_patch(fout, patch, tier, mult=10)
        left -= 5

plt.savefig("patches.png")
