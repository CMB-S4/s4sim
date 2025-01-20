#  Generate a list of observing targets to match SAT depth

import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np


fname_out = "patches.txt"

"""
cov = hp.read_map(
    "/global/cfs/cdirs/cmbs4/chile_optimization/simulations/phase2/noise_depth/sun90_f155_90years_cov.fits"
)
rhits = np.zeros_like(cov)
good = cov != 0
rhits[good] = 1 / cov[good]
"""
rhits = hp.read_map("../../sims/outputs/sun90max/f085/season/mapmaker_hits.fits", dtype=float)
rhits += hp.read_map("../../sims/outputs/sun90max/f085/break/mapmaker_hits.fits")
rhits /= np.amax(rhits)
rhits = rhits**2  # Use the square of the hits as weight

nside = hp.get_nside(rhits)
npix = rhits.size
pix = np.arange(npix)
lons, lats = hp.pix2ang(nside, pix, lonlat=True)

throw = 20.0
scantime = 20.0

dlon = 5
dlat = 1
r = np.radians(5)
lines = []
for lat in range(-90, 90, dlat):
    for lon in range(0, 360, dlon):
        vec = hp.ang2vec(lon, lat, lonlat=True)
        pix = hp.query_disc(nside, vec, radius=r)
        hitmean = np.mean(rhits[pix])
        if hitmean < 1e-2:
            continue
        # Only high priority targets can be targets with the edge of the focalplane
        if hitmean > 0.5:
            radius = 3.0
        else:
            radius = 0.1
        priority = 1 / hitmean
        if lon >= 100 and lon <= 190:
            # Lower priority for the Northern tiles
            priority *= 1000
        name = f"DEC{lat:+04}_RA{lon:+04}"
        lines.append(
            f"{name},MAX-DEPTH,{priority:.3e},{lon:.3f},{lat:.3f},{radius},{throw},{scantime}"
        )

with open(fname_out, "w") as f:
    for line in lines:
        f.write("--patch\n")
        f.write(line + "\n")
print(f"Wrote {len(lines)} targets to {fname_out}")
