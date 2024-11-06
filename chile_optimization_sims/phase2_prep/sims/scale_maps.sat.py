# This script scales the simulated noise covariances with efficiency and
# survey factors that were not included in the simulation

import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np


for isplit in range(8, 10):
    split = f"split{isplit}"

    # Total target efficiency factor

    f_total = {}

    # fsky for measurement requirement

    fsky = 0.03

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

    # CHSAT f_total from
    # https://docs.google.com/spreadsheets/d/116Xa1vHrIwO6xTLsZalXo-QK7aKRQJnTE5LhkSl9eig/edit?gid=950965928#gid=950965928
    # on 10/07/2024

    f_total["sat"] = {
        "f030" : 0.17,
        "f040" : 0.17,
        "f085" : 0.17,
        "f090" : 0.17,
        "f095" : 0.17,
        "f145" : 0.16,
        "f150" : 0.16,
        "f155" : 0.16,
        "f220" : 0.10,
        "f280" : 0.07,
    }
    f_field["sat"] = 0.54

    # f_weather from the simulation logs on 11/06/2024

    if isplit in [0, 3, 6, 7, 8]:
        cut_3mm = 1 - 955 / 5452
        cut_2mm = 1 - 1636 / 5452
        f_weather_sim["sat"] = {
            "f030" : cut_3mm,
            "f040" : cut_3mm,
            "f085" : cut_3mm,
            "f090" : cut_3mm,
            "f095" : cut_3mm,
            "f145" : cut_3mm,
            "f150" : cut_3mm,
            "f155" : cut_3mm,
            "f220" : cut_2mm,
            "f280" : cut_2mm,
        }
    elif isplit in [1, 4]:
        cut_3mm = 1 - 1234 / 7078
        cut_2mm = 1 - 2100 / 7078
        f_weather_sim["sat"] = {
            "f030" : cut_3mm,
            "f040" : cut_3mm,
            "f085" : cut_3mm,
            "f090" : cut_3mm,
            "f095" : cut_3mm,
            "f145" : cut_3mm,
            "f150" : cut_3mm,
            "f155" : cut_3mm,
            "f220" : cut_2mm,
            "f280" : cut_2mm,
        }
    elif isplit in [2, 5, 9]:
        cut_3mm = 1 - 1452 / 8087
        cut_2mm = 1 - 2473 / 8087
        f_weather_sim["sat"] = {
            "f030" : cut_3mm,
            "f040" : cut_3mm,
            "f085" : cut_3mm,
            "f090" : cut_3mm,
            "f095" : cut_3mm,
            "f145" : cut_3mm,
            "f150" : cut_3mm,
            "f155" : cut_3mm,
            "f220" : cut_2mm,
            "f280" : cut_2mm,
        }
    else:
        raise RuntimeError(f"Unknown split: {isplit}")


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

    for flavor in "sat",:
        nrow, ncol = 2, 4
        fig = plt.figure(figsize=[4 * ncol, 4 * nrow])
        fig.suptitle(f"{flavor} {split}, {nsat * 3} tubes, fsky = {fsky}")
        iplot = 0
        for band in f_total[flavor].keys():
            fname_in = f"outputs/{flavor}/{split}/{band}/mapmaker_cov.fits"
            if not os.path.isfile(fname_in):
                continue
            print(f"Loading {fname_in}")
            cov = hp.read_map(fname_in, None)

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
            scale *= f_weather_sim[flavor][band]
            scale *= f_season
            scale *= f_field[flavor]
            # Now apply f_total
            scale /= f_total[flavor][band]
            # Scale
            cov *= scale

            # Save the scaled covariance

            outdir = f"scaled_outputs/sat/{split}"
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
            depth = depth_I * np.sqrt(2)  # depth_P
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
        fname_plot = f"plots/{flavor}_{split}.png"
        fig.savefig(fname_plot)
        print(f"Plot saved in {fname_plot}")
