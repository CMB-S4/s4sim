import os
import pickle
import sys

import astropy.units as u
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import c, h, k

from toast.pixels_io_healpix import read_healpix


nrow, ncol = 1, 1
fig = plt.figure(figsize=[6 * ncol, 4 * nrow])
ax1 = fig.add_subplot(nrow, ncol, 1)

total_depth = {}
for band in "f090", "f150":
    fname_in = f"depths_splat_{band}.pck"
    if os.path.isfile(fname_in):
        with open(fname_in, "rb") as f:
            depths = pickle.load(f)
    else:
        depths = {}

    x, y = [], []
    for key, value in depths.items():
        x.append(key)
        y.append(value)
        if key in total_depth:
            depth1 = total_depth[key]
            depth2 = value
            total_depth[key] = 1 / np.sqrt(1 / depth1**2 + 1 / depth2**2)
        else:
            total_depth[key] = value
    y = np.array(y)
    doy = np.arange(y.size) + 1
    ax1.plot(doy, y, '.', label=f"{band}")

x, y = [], []
for key, value in total_depth.items():
    x.append(key)
    y.append(value)
y = np.array(y)
doy = np.arange(y.size) + 1
ax1.plot(doy, y, '.', label=f"Combined")

ax1.set_xlabel("DOY")
ax1.set_ylabel("Depth [mJy]")
ax1.axhline(3.0, color="k", linestyle="--", label="MR 4.2")
ax1.set_title("Daily SPLAT depth")
ax1.legend(bbox_to_anchor=(1.0, 1.0))
ax1.set_ylim([0, 5])
ax1.set_xlim([0, 365])

fig.tight_layout()
plt.savefig("splat_depth.png")
plt.savefig("splat_depth.pdf")
