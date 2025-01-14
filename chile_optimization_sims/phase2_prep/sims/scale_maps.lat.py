# This script scales the simulated noise covariances with efficiency and
# survey factors that were not included in the simulation

import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np


for isplit in [1, 2]:
    split = f"split{isplit}"
    flavors = {
        "split1" : ["lat_wide", "lat_delensing"],
        "split2" : ["lat_delensing"],
    }[split]

    # Total target efficiency factor

    f_total = {}

    # fsky for measurement requirement

    fskies = {}

    # Factors not included in f_total

    yield_ = 0.8
    f_weight = 1.0

    # Efficiency factors already applied in the simulation must
    # not be double-counted in f_total

    f_season = 0.75
    f_field = {}

    # Weather cut is approximately included in the sim but we want to remove
    # it and replace with the factor included in f_total

    f_weather_sim = {}

    # LAT wide f_total from
    # https://docs.google.com/spreadsheets/d/116Xa1vHrIwO6xTLsZalXo-QK7aKRQJnTE5LhkSl9eig/edit?gid=226211982#gid=226211982
    # on 09/27/2024

    f_total["lat_wide"] = {
        "f020" : 0.29,
        "f030" : 0.29,
        "f040" : 0.29,
        "f090" : 0.29,
        "f150" : 0.29,
        "f220" : 0.25,
        "f280" : 0.18,
    }
    f_field["lat_wide"] = 0.908
    fskies["lat_wide"] = 0.6

    # LAT delensing f_total from
    # https://docs.google.com/spreadsheets/d/116Xa1vHrIwO6xTLsZalXo-QK7aKRQJnTE5LhkSl9eig/edit?gid=226211982#gid=226211982
    # on 09/27/2024

    f_total["lat_delensing"] = {
        "f020" : 0.30,
        "f030" : 0.30,
        "f040" : 0.30,
        "f090" : 0.30,
        "f150" : 0.30,
        "f220" : 0.26,
        "f280" : 0.18,
    }
    f_field["lat_delensing"] = 0.972
    fskies["lat_delensing"] = 0.03

    # f_weather from the simulation logs on 11/06/2024

    cut_3mm = 1 - 1218 / 6629
    cut_2mm = 1 - 2067 / 6629
    f_weather_sim["lat_wide"] = {
        "f020" : cut_3mm,
        "f030" : cut_3mm,
        "f040" : cut_3mm,
        "f090" : cut_3mm,
        "f150" : cut_3mm,
        "f220" : cut_2mm,
        "f280" : cut_2mm,
    }
    if split == "split1":
        cut_3mm = 1 - 1535 / 8494
        cut_2mm = 1 - 2640 / 8494
        f_weather_sim["lat_delensing"] = {
            "f020" : cut_3mm,
            "f030" : cut_3mm,
            "f040" : cut_3mm,
            "f090" : cut_3mm,
            "f150" : cut_3mm,
            "f220" : cut_2mm,
            "f280" : cut_2mm,
        }
    elif split == "split2":
        cut_3mm = 1 - 1549 / 8601
        cut_2mm = 1 - 2645 / 8601
        f_weather_sim["lat_delensing"] = {
            "f020" : cut_3mm,
            "f030" : cut_3mm,
            "f040" : cut_3mm,
            "f090" : cut_3mm,
            "f150" : cut_3mm,
            "f220" : cut_2mm,
            "f280" : cut_2mm,
        }
    else:
        raise RuntimeError(f"Unknown split: {split}")

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

    # Number of telescope years is a parameter

    n_years = {
        "lat_wide" : [14],
        "lat_delensing" : [16, 26, 36],
    }


    # Loop over all covariance matrices

    for flavor in flavors:
        for n_year in n_years[flavor]:
            nrow, ncol = 2, 4
            fig = plt.figure(figsize=[4 * ncol, 4 * nrow])
            fig.suptitle(f"{flavor} {split} {n_year} years")
            iplot = 0
            for band in f_total[flavor].keys():
                fname_in = f"outputs/{flavor}/{split}/{band}/mapmaker_cov.fits"
                print(f"Loading {fname_in}")
                cov = hp.read_map(fname_in, None)

                # Scale the white noise covariance

                scale = 1.
                # Compensate for focalplane decimation
                scale /= thinfp[band]
                # Account for full mission length
                scale /= n_year
                # Yield and f_weight are not in f_total
                scale /= yield_
                scale /= f_weight
                # Eliminate factors from simulation (included in f_total)
                scale *= f_weather_sim[flavor][band]
                scale *= f_season
                scale *= f_field[flavor]
                # Now apply f_total
                scale /= f_total[flavor][band]
                # Scale
                cov *= scale

                # Save the scaled covariance

                outdir = f"scaled_outputs/lat/{split}"
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
                #
                sorted_depth = depth.copy()
                sorted_depth[sorted_depth == 0] = 1e10
                sorted_depth = np.sort(sorted_depth)
                fsky = fskies[flavor]
                lim = int(depth.size * fsky)
                mean_depth = np.mean(sorted_depth[:lim])
                #
                depth[depth == 0] = hp.UNSEEN
                hp.mollview(
                    depth,
                    min=vmin,
                    max=vmax,
                    title=f"{band}, mean = {mean_depth:.2f} (fsky={fsky})",
                    sub=[nrow, ncol, iplot],
                    cmap="inferno",
                    unit="$\mu$K.arcmin",
                    xsize=1600,
                    format="%.3f"
                )

            # Save plot

            os.makedirs("plots", exist_ok=True)
            fname_plot = f"plots/{flavor}_{split}_{n_year}years.png"
            fig.savefig(fname_plot)
            print(f"Plot saved in {fname_plot}")
