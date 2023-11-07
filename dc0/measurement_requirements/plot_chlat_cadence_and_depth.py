from datetime import date, timedelta
import os
import pickle
import sys

import astropy.units as u
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import c, h, k

from toast.pixels_io_healpix import read_healpix


nrow, ncol = 1, 2
fig = plt.figure(figsize=[6 * ncol, 4 * nrow])
ax1 = fig.add_subplot(nrow, ncol, 1)
ax2 = fig.add_subplot(nrow, ncol, 2)

total_depth = {}
for band in "f090", "f150":
    fname_in = f"cadence_and_depth_chlat_{band}.pck"
    with open(fname_in, "rb") as f:
        common, depths = pickle.load(f)

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
ax1.axhline(10.0, color="k", linestyle="--", label="MR 4.2")
ax1.axvline(90.0, color="k", linestyle="-", label="Start of season")
ax1.set_title("Daily CHLAT depth")
#ax1.legend(bbox_to_anchor=(1.0, 1.0))
ax1.legend(loc="best")
ax1.set_ylim([0, 20])
ax1.set_xlim([0, 365])

x, y = [], []
for key, value in common.items():
    x.append(key)
    y.append(value)
y = np.array(y)
doy = np.arange(y.size) + 1
ax2.plot(doy, y, '.', label="Last 5 days")
ax2.set_xlabel("DOY")
ax2.set_ylabel("fsky")
ax2.axhline(0.25, color="k", linestyle="--", label="MR 4.1")
ax2.set_title("5-day common sky fraction")
ax2.legend(loc="best")

fig.tight_layout()
plt.savefig("chlat_cadence_and_depth.png")
plt.savefig("chlat_cadence_and_depth.pdf")
