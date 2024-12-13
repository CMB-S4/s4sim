# This script scales the simulated noise covariances with efficiency and
# survey factors that were not included in the simulation

import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np


# Total target efficiency factor

f_total = {}

# fsky for measurement requirement

fsky = 0.03

# Factors not included in f_total

yield_ = 0.8
f_weight = 1.0

# Efficiency factors already applied in the simulation must
# not be double-counted in f_total

f_season = {
    "season" : 0.75,
    "break" : 0.25,  # Cannot use 0.05 or will double-correct
}
f_field = {
    "sun90" : {
        "season" : 0.57,
        "break" : 0.53,
    },
    "sun45" : {
        "season" : 0.83,
        "break" : 0.76,
    },
}

# Weather cut is approximately included in the sim but we want to remove
# it and replace with the factor included in f_total

f_weather_sim = {}

# CHSAT f_total from
# https://docs.google.com/spreadsheets/d/1QtA3dp0GGRWWO1WEGKVnvEYFx-MhwpxRN3enGbCDOi8/edit?pli=1&gid=774702783#gid=774702783
# on 12/11/2024

f_sensitivity = yield_ * f_weight
f_total = {
    "sun90" : {
        "season" : {
            "f030" : 0.140 / f_sensitivity,
            "f040" : 0.140 / f_sensitivity,
            "f085" : 0.140 / f_sensitivity,
            "f090" : 0.140 / f_sensitivity,
            "f095" : 0.140 / f_sensitivity,
            "f145" : 0.140 / f_sensitivity,
            "f150" : 0.140 / f_sensitivity,
            "f155" : 0.140 / f_sensitivity,
            "f220" : 0.086 / f_sensitivity,
            "f280" : 0.060 / f_sensitivity,
        },
        "break" : {
            "f030" : 0.011 / f_sensitivity,
            "f040" : 0.011 / f_sensitivity,
            "f085" : 0.011 / f_sensitivity,
            "f090" : 0.011 / f_sensitivity,
            "f095" : 0.011 / f_sensitivity,
            "f145" : 0.010 / f_sensitivity,
            "f150" : 0.010 / f_sensitivity,
            "f155" : 0.010 / f_sensitivity,
            "f220" : 0.006 / f_sensitivity,
            "f280" : 0.004 / f_sensitivity,
        }
    },
    "sun45" : {
        "season" : {
            "f030" : 0.203 / f_sensitivity,
            "f040" : 0.203 / f_sensitivity,
            "f085" : 0.203 / f_sensitivity,
            "f090" : 0.203 / f_sensitivity,
            "f095" : 0.203 / f_sensitivity,
            "f145" : 0.200 / f_sensitivity,
            "f150" : 0.200 / f_sensitivity,
            "f155" : 0.200 / f_sensitivity,
            "f220" : 0.125 / f_sensitivity,
            "f280" : 0.087 / f_sensitivity,
        },
        "break" : {
            "f030" : 0.015 / f_sensitivity,
            "f040" : 0.015 / f_sensitivity,
            "f085" : 0.015 / f_sensitivity,
            "f090" : 0.015 / f_sensitivity,
            "f095" : 0.015 / f_sensitivity,
            "f145" : 0.015 / f_sensitivity,
            "f150" : 0.015 / f_sensitivity,
            "f155" : 0.015 / f_sensitivity,
            "f220" : 0.008 / f_sensitivity,
            "f280" : 0.006 / f_sensitivity,
        }
    },
}

# f_weather from prune_schedule.py

cut_3mm_season = 6631 / 8099
cut_2mm_season = 5593 / 8099
cut_3mm_break = 498 / 2495
cut_2mm_break = 378 / 2495
f_weather_sim["sun45"] = {
    "season" : {
        "f030" : cut_3mm_season,
        "f040" : cut_3mm_season,
        "f085" : cut_3mm_season,
        "f090" : cut_3mm_season,
        "f095" : cut_3mm_season,
        "f145" : cut_3mm_season,
        "f150" : cut_3mm_season,
        "f155" : cut_3mm_season,
        "f220" : cut_2mm_season,
        "f280" : cut_2mm_season,
    },
    "break" : {
        "f030" : cut_3mm_break,
        "f040" : cut_3mm_break,
        "f085" : cut_3mm_break,
        "f090" : cut_3mm_break,
        "f095" : cut_3mm_break,
        "f145" : cut_3mm_break,
        "f150" : cut_3mm_break,
        "f155" : cut_3mm_break,
        "f220" : cut_2mm_break,
        "f280" : cut_2mm_break,
    },
}

cut_3mm_season = 4429 / 5476
cut_2mm_season = 3840 / 5476
cut_3mm_break = 412 / 1612
cut_2mm_break = 267 / 1612
f_weather_sim["sun90"] = {
    "season" : {
        "f030" : cut_3mm_season,
        "f040" : cut_3mm_season,
        "f085" : cut_3mm_season,
        "f090" : cut_3mm_season,
        "f095" : cut_3mm_season,
        "f145" : cut_3mm_season,
        "f150" : cut_3mm_season,
        "f155" : cut_3mm_season,
        "f220" : cut_2mm_season,
        "f280" : cut_2mm_season,
    },
    "break" : {
        "f030" : cut_3mm_break,
        "f040" : cut_3mm_break,
        "f085" : cut_3mm_break,
        "f090" : cut_3mm_break,
        "f095" : cut_3mm_break,
        "f145" : cut_3mm_break,
        "f150" : cut_3mm_break,
        "f155" : cut_3mm_break,
        "f220" : cut_2mm_break,
        "f280" : cut_2mm_break,
    },
}

# Focalplane decimation factors must be compensated for

thinfp = 1

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
    "f090" : nyear * ntube_mf * 2,
    "f095" : nyear * ntube_mf,
    "f145" : nyear * ntube_mf,
    "f150" : nyear * ntube_mf * 2,
    "f155" : nyear * ntube_mf,
    "f220" : nyear * ntube_hf,
    "f280" : nyear * ntube_hf,
}


# Loop over all covariance matrices

for flavor in "sun90", "sun45":
    nrow, ncol = 2, 5
    fig = plt.figure(figsize=[4 * ncol, 4 * nrow])
    fig.suptitle(f"{flavor}, {nsat * 3} tubes, fsky = {fsky}")
    iplot = 0
    for band in f_total[flavor]["season"].keys():
        covs = []
        for period in "season", "break":
            fname_in = f"outputs/{flavor}/{band}/{period}/mapmaker_cov.fits"
            print(f"Loading {fname_in}")
            cov = hp.read_map(fname_in, [0, 3, 5])  # II, QQ, UU

            # Scale the white noise covariance

            n_year = n_years[band]

            scale = 1.
            # Compensate for focalplane decimation
            scale /= thinfp
            # Account for full mission length
            scale /= n_year
            # Yield and f_weight are not in f_total
            scale /= yield_
            scale /= f_weight
            # Eliminate factors from simulation (included in f_total)
            scale *= f_weather_sim[flavor][period][band]
            scale *= f_season[period]
            scale *= f_field[flavor][period]
            # Now apply f_total
            scale /= f_total[flavor][period][band]
            # Scale
            cov *= scale
            covs.append(cov)

        # Combine the covariances
        invcov_sum = None
        for cov in covs:
            invcov = np.zeros_like(cov)
            invcov[cov != 0] = 1 / cov[cov != 0]
            if invcov_sum is None:
                invcov_sum = invcov
            else:
                invcov_sum += invcov
        cov = np.zeros_like(invcov_sum)
        good = invcov_sum != 0
        cov[good] = 1 / invcov_sum[good]

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
        depth_Q = np.sqrt(cov[1] * pix_area) * 1e6  # uK.arcmin
        depth_U = np.sqrt(cov[2] * pix_area) * 1e6  # uK.arcmin
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
        #
        sorted_depth = depth.copy()
        sorted_depth[sorted_depth == 0] = 1e10
        sorted_depth = np.sort(sorted_depth)
        lim = int(depth.size * fsky)
        mean_depth = np.mean(sorted_depth[:lim])
        #
        depth[depth == 0] = hp.UNSEEN
        hp.mollview(
            depth,
            min=vmin,
            max=vmax,
            title=f"{band}, {n_year} tube years, mean = {mean_depth:.2f}",
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
