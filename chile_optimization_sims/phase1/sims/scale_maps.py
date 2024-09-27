# This script scales the simulated noise covariances with efficiency and
# survey factors that were not included in the simulation

import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np


# Total target efficiency factor

f_total = {}

# Factors not included in f_total

yield_ = 0.8
f_weight = 0.8
n_year = 10

# Efficiency factors already applied in the simulation must
# not be double-counted in f_total

f_season = 0.75
f_field = {}
f_scanset = {}

# Weather cut is approximately include in the sim but we want to remove
# it and replace with the factor included in f_total
f_weather_sim = {}

# LAT wide f_total from
# https://docs.google.com/spreadsheets/d/116Xa1vHrIwO6xTLsZalXo-QK7aKRQJnTE5LhkSl9eig/edit?gid=226211982#gid=226211982
# on 09/27/2024

f_total["lat_wide"] = {
    "f020" : 0.31,
    "f030" : 0.31,
    "f040" : 0.31,
    "f090" : 0.31,
    "f150" : 0.31,
    "f220" : 0.27,
    "f280" : 0.19,
}
f_field["lat_wide"] = 0.979
f_scanset["lat_wide"] = 0.92

# LAT delensing f_total from
# https://docs.google.com/spreadsheets/d/116Xa1vHrIwO6xTLsZalXo-QK7aKRQJnTE5LhkSl9eig/edit?gid=226211982#gid=226211982
# on 09/27/2024

f_total["lat_delensing"] = {
    "f020" : 0.30,
    "f030" : 0.30,
    "f040" : 0.30,
    "f090" : 0.30,
    "f150" : 0.30,
    "f220" : 0.27,
    "f280" : 0.18,
}
f_field["lat_delensing"] = 0.998
f_scanset["lat_delensing"] = 0.92

# LAT delensing_core is a variant of LAT delensing

f_total["lat_delensing_core"] = f_total["lat_delensing"]
f_field["lat_delensing_core"] = f_field["lat_delensing"]
f_scanset["lat_delensing_core"] = f_scanset["lat_delensing"]

# f_weather from the simulation logs on 09/27/2024

f_weather_sim["lat_wide"] = {
    "f020" : 1 - 1163 / 6214,
    "f030" : 1 - 1163 / 6214,
    "f040" : 1 - 1163 / 6214,
    "f090" : 1 - 1163 / 6214,
    "f150" : 1 - 1163 / 6214,
    "f220" : 1 - 1968 / 6214,
    "f280" : 1 - 1968 / 6214,
}
f_weather_sim["lat_delensing"] = {
    "f020" : 1 - 1533 / 8366,
    "f030" : 1 - 1533 / 8366,
    "f040" : 1 - 1533 / 8366,
    "f090" : 1 - 1533 / 8366,
    "f150" : 1 - 1533 / 8366,
    "f220" : 1 - 2621 / 8366,
    "f280" : 1 - 2621 / 8366,
}
f_weather_sim["lat_delensing_core"] = {
    "f020" : 1 - 1538 / 8438,
    "f030" : 1 - 1538 / 8438,
    "f040" : 1 - 1538 / 8438,
    "f090" : 1 - 1538 / 8438,
    "f150" : 1 - 1538 / 8438,
    "f220" : 1 - 2639 / 8438,
    "f280" : 1 - 2639 / 8438,
}

# Focalplane decimation factors must be compensated for

thinfp = {
    "f020" : 1,
    "f030" : 1,
    "f040" : 1,
    "f090" : 16,
    "f150" : 16,
    "f220" : 16,
    "f280" : 16,
}


# Loop over all covariance matrices

for flavor in "lat_wide", "lat_delensing", "lat_delensing_core":
    nrow, ncol = 2, 4
    fig = plt.figure(figsize=[4 * ncol, 4 * nrow])
    fig.suptitle(flavor)
    iplot = 0
    telescope = flavor[:3]
    for band in f_total[flavor].keys():
        fname_in = f"outputs/{telescope}/{flavor}/{band}/mapmaker_cov.fits"
        print(f"Loading {fname_in}")
        cov = hp.read_map(fname_in, None)

        # Scale the white noise covariance

        cov /= thinfp[band]
        cov /= yield_
        cov /= f_weight
        cov /= n_year
        cov /= (
            f_total[flavor][band]
            / f_season
            / f_field[flavor]
            / f_scanset[flavor]
            * f_weather_sim[flavor][band]
        )

        # Save the scaled covariance

        outdir = f"scaled_outputs"
        os.makedirs(outdir, exist_ok=True)
        fname_out = f"{outdir}/{flavor}_{band}_cov.fits"
        print(f"Writing {fname_out}")
        hp.write_map(fname_out, cov, dtype=np.float32, coord="C", overwrite=True)

        # Derive depth from the covariance and save

        nside = hp.get_nside(cov)
        pix_area = hp.nside2pixarea(nside, degrees=True) * 3600  # arcmin^2
        depth_I = np.sqrt(cov[0] * pix_area) * 1e6  # uK.arcmin
        depth_Q = np.sqrt(cov[3] * pix_area) * 1e6  # uK.arcmin
        depth_U = np.sqrt(cov[5] * pix_area) * 1e6  # uK.arcmin
        fname_out = f"{outdir}/{flavor}_{band}_depth.fits"
        print(f"Writing {fname_out}")
        hp.write_map(fname_out, [depth_I, depth_Q, depth_U], dtype=np.float32, coord="C", overwrite=True)

        # Plot depth

        iplot += 1
        depth = depth_I
        vmin = np.amin(depth[depth != 0])
        vmax = 2 * vmin
        depth[depth == 0] = hp.UNSEEN
        hp.mollview(depth, min=vmin, max=vmax, title=band, sub=[nrow, ncol, iplot], cmap="inferno")

    # Save plot

    os.makedirs("plots", exist_ok=True)
    fname_plot = f"plots/{flavor}.png"
    fig.savefig(fname_plot)
    print(f"Plot saved in {fname_plot}")
