
# This script scales the simulated noise covariances with efficiency and
# survey factors that were not included in the simulation

import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np

# Factors copied from
# https://docs.google.com/spreadsheets/d/1n2NyRSKN9OZRtLJp6FTTG66upSJUAcfMDarWwU9IYb0/edit?pli=1&gid=859341864#gid=859341864
# on 2025/01/17


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
    "sun90max" : {
        "season" : 0.37,
        "break" : 0.40,
    },
    "sun45_supplement" : {
        "season" : 0.56,
        "break" : 0.64,
    },
}

# Weather cut is approximately included in the sim but we want to remove
# it and replace with the factor included in f_total

f_weather_sim = {}

f_sensitivity = yield_ * f_weight
f_total = {
    "sun90max" : {
        "season" : {
            "f030" : 0.091 / f_sensitivity,
            "f040" : 0.091 / f_sensitivity,
            "f085" : 0.091 / f_sensitivity,
            "f090" : 0.091 / f_sensitivity,
            "f095" : 0.091 / f_sensitivity,
            "f145" : 0.090 / f_sensitivity,
            "f150" : 0.090 / f_sensitivity,
            "f155" : 0.090 / f_sensitivity,
            "f220" : 0.056 / f_sensitivity,
            "f280" : 0.039 / f_sensitivity,
        },
        "break" : {
            "f030" : 0.008 / f_sensitivity,
            "f040" : 0.008 / f_sensitivity,
            "f085" : 0.008 / f_sensitivity,
            "f090" : 0.008 / f_sensitivity,
            "f095" : 0.008 / f_sensitivity,
            "f145" : 0.008 / f_sensitivity,
            "f150" : 0.008 / f_sensitivity,
            "f155" : 0.008 / f_sensitivity,
            "f220" : 0.008 / f_sensitivity,
            "f280" : 0.008 / f_sensitivity,
        }
    },
    "sun45_supplement" : {
        "season" : {
            "f030" : 0.138 / f_sensitivity,
            "f040" : 0.138 / f_sensitivity,
            "f085" : 0.138 / f_sensitivity,
            "f090" : 0.138 / f_sensitivity,
            "f095" : 0.138 / f_sensitivity,
            "f145" : 0.136 / f_sensitivity,
            "f150" : 0.136 / f_sensitivity,
            "f155" : 0.136 / f_sensitivity,
            "f220" : 0.085 / f_sensitivity,
            "f280" : 0.059 / f_sensitivity,
        },
        "break" : {
            "f030" : 0.013 / f_sensitivity,
            "f040" : 0.013 / f_sensitivity,
            "f085" : 0.013 / f_sensitivity,
            "f090" : 0.013 / f_sensitivity,
            "f095" : 0.013 / f_sensitivity,
            "f145" : 0.013 / f_sensitivity,
            "f150" : 0.013 / f_sensitivity,
            "f155" : 0.013 / f_sensitivity,
            "f220" : 0.006 / f_sensitivity,
            "f280" : 0.004 / f_sensitivity,
        }
    },
}

# f_weather from prune_schedule.py

cut_3mm_season = 3006 / 3703
cut_2mm_season = 2504 / 3703
cut_3mm_break = 296 / 1385
cut_2mm_break = 209 / 1385
f_weather_sim["sun45_supplement"] = {
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

cut_3mm_season = 1977 / 2454
cut_2mm_season = 1661 / 2454
cut_3mm_break = 233 / 854
cut_2mm_break = 140 / 854
f_weather_sim["sun90max"] = {
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

for flavor in "sun90max", "sun45_supplement":
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
    fname_plot = f"plots/{flavor}.png"
    fig.savefig(fname_plot)
    print(f"Plot saved in {fname_plot}")
