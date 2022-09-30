from glob import glob
import os
import sys

from astropy.io.misc.hdf5 import read_table_hdf5
import matplotlib.pyplot as plt
import numpy as np


pwv = np.linspace(0, 3, 100)
el = np.linspace(40, 60, 100)

nrow, ncol = 2, 1
fig = plt.figure(figsize=[8 * ncol, 4 * nrow])
ax1 = fig.add_subplot(nrow, ncol, 1)
ax2 = fig.add_subplot(nrow, ncol, 2)
ax1.set_title("NET ratio vs. elevation")
ax2.set_title("NET ratio vs. PWV")
ax1.set_xlabel("Elevation [deg]")
ax2.set_xlabel("PWV [mm]")

fmt = "--"
for band in "f030", "f040", "f085", "f095", "f145", "f155", "f220", "f280":
    if band == "f085":
        fmt = "-"
    elif band == "f220":
        fmt = ":"
    table_pole = read_table_hdf5(glob(f"focalplane_SAT*_SAT_{band}.h5")[0])
    table_chile = read_table_hdf5(glob(f"focalplane_CHSAT*_CHSAT_{band}.h5")[0])

    net_pole = table_pole["psd_net"][0]
    pwv_a0_pole = table_pole["pwv_noise_a0"][0]
    pwv_a1_pole = table_pole["pwv_noise_a1"][0]
    pwv_a2_pole = table_pole["pwv_noise_a2"][0]
    elevation_a_pole = table_pole["elevation_noise_a"][0]
    elevation_c_pole = table_pole["elevation_noise_c"][0]

    net_chile = table_chile["psd_net"][0]
    pwv_a0_chile = table_chile["pwv_noise_a0"][0]
    pwv_a1_chile = table_chile["pwv_noise_a1"][0]
    pwv_a2_chile = table_chile["pwv_noise_a2"][0]
    elevation_a_chile = table_chile["elevation_noise_a"][0]
    elevation_c_chile = table_chile["elevation_noise_c"][0]

    pwv_factor_pole = pwv_a0_pole + pwv_a1_pole * pwv + pwv_a2_pole * pwv**2
    elevation_factor_pole = elevation_a_pole / np.sin(np.radians(el)) + elevation_c_pole

    pwv_factor_chile = pwv_a0_chile + pwv_a1_chile * pwv + pwv_a2_chile * pwv**2
    elevation_factor_chile = elevation_a_chile / np.sin(np.radians(el)) + elevation_c_chile

    ax1.plot(el, net_chile * elevation_factor_chile / net_pole / elevation_factor_pole, fmt, label=band)
    ax2.plot(pwv, net_chile * pwv_factor_chile / net_pole / pwv_factor_pole, fmt, label=band)

plt.legend(loc="best")
plt.tight_layout()
plt.savefig("net_comparison.png")
