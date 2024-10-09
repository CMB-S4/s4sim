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
f_weight = {
    "f030" : 0.8,
    "f040" : 0.8,
    "f085" : 0.8,
    "f095" : 0.8,
    "f145" : 0.65,
    "f155" : 0.65,
    "f220" : 0.65,
    "f280" : 0.65,
}

# Efficiency factors already applied in the simulation must
# not be double-counted in f_total

f_season = 0.75
f_field = {}

# Weather cut is approximately included in the sim but we want to remove
# it and replace with the factor included in f_total

f_weather_sim = {}

# CHSAT f_total from
# https://docs.google.com/spreadsheets/d/116Xa1vHrIwO6xTLsZalXo-QK7aKRQJnTE5LhkSl9eig/edit?gid=950965928#gid=950965928
# on 10/07/2024

f_total["sat"] = {
    "f030" : 0.17,
    "f040" : 0.17,
    "f085" : 0.17,
    "f095" : 0.17,
    "f145" : 0.16,
    "f155" : 0.16,
    "f220" : 0.10,
    "f280" : 0.07,
}
f_field["sat"] = 0.54

# f_weather from the simulation logs on 10/07/2024

cut_3mm = 1 - 903 / 5157
cut_2mm = 1 - 1525 / 5157
f_weather_sim["sat"] = {
    "f030" : cut_3mm,
    "f040" : cut_3mm,
    "f085" : cut_3mm,
    "f095" : cut_3mm,
    "f145" : cut_3mm,
    "f155" : cut_3mm,
    "f220" : cut_2mm,
    "f280" : cut_2mm,
}

# Focalplane decimation factors must be compensated for

thinfp = {
    "f030" : 1,
    "f040" : 1,
    "f085" : 1,
    "f095" : 1,
    "f145" : 1,
    "f155" : 1,
    "f220" : 1,
    "f280" : 1,
}

# Number of telescope years is a parameter

nsat = 9
ntube_lf = nsat // 3
ntube_mf = nsat
ntube_hf = nsat // 3 * 2
nyear = 10
n_years = {
    "f030" : nyear * ntube_lf,
    "f040" : nyear * ntube_lf,
    "f085" : nyear * ntube_mf,
    "f095" : nyear * ntube_mf,
    "f145" : nyear * ntube_mf,
    "f155" : nyear * ntube_mf,
    "f220" : nyear * ntube_hf,
    "f280" : nyear * ntube_hf,
}


# Loop over all covariance matrices

for flavor in "sat",:
    nrow, ncol = 2, 4
    fig = plt.figure(figsize=[4 * ncol, 4 * nrow])
    fig.suptitle(f"{flavor}, {nsat * 3} tubes")
    iplot = 0
    for band in f_total[flavor].keys():
        fname_in = f"outputs/{flavor}/{band}/mapmaker_cov.fits"
        print(f"Loading {fname_in}")
        cov = hp.read_map(fname_in, None)

        # Scale the white noise covariance

        n_year = n_years[band]

        scale = 1.
        # Compensate for focalplane decimation
        scale /= thinfp[band]
        # Account for full mission length
        scale /= n_year
        # Yield and f_weight are not in f_total
        scale /= yield_
        scale /= f_weight[band]
        # Eliminate factors from simulation (included in f_total)
        scale *= f_weather_sim[flavor][band]
        scale *= f_season
        scale *= f_field[flavor]
        # Now apply f_total
        scale /= f_total[flavor][band]
        # Scale
        cov *= scale

        # Save the scaled covariance

        outdir = f"scaled_outputs"
        os.makedirs(outdir, exist_ok=True)
        fname_out = f"{outdir}/{flavor}_{band}_{n_year}years_cov.fits"
        print(f"Writing {fname_out}")
        hp.write_map(fname_out, cov, dtype=np.float32, coord="C", overwrite=True)

        # Derive depth from the covariance and save

        nside = hp.get_nside(cov)
        pix_area = hp.nside2pixarea(nside, degrees=True) * 3600  # arcmin^2
        depth_I = np.sqrt(cov[0] * pix_area) * 1e6  # uK.arcmin
        depth_Q = np.sqrt(cov[3] * pix_area) * 1e6  # uK.arcmin
        depth_U = np.sqrt(cov[5] * pix_area) * 1e6  # uK.arcmin
        fname_out = f"{outdir}/{flavor}_{band}_{n_year}years_depth.fits"
        print(f"Writing {fname_out}")
        hp.write_map(
            fname_out,
            [depth_I, depth_Q, depth_U],
            dtype=np.float32,
            coord="C",
            overwrite=True,
        )

        # Plot depth

        iplot += 1
        depth = depth_I
        vmin = np.amin(depth[depth != 0])
        vmax = 2 * vmin
        depth[depth == 0] = hp.UNSEEN
        hp.mollview(
            depth,
            min=vmin,
            max=vmax,
            title=f"{band}, {n_year} tube years",
            sub=[nrow, ncol, iplot],
            cmap="inferno",
            unit="$\mu$K.arcmin",
            xsize=1600,
            format="%.3f"
        )

    # Save plot

    os.makedirs("plots", exist_ok=True)
    fname_plot = f"plots/{flavor}.png"
    fig.savefig(fname_plot)
    print(f"Plot saved in {fname_plot}")
