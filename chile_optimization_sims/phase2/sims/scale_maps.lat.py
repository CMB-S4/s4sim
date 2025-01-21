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

fskies = {}

# Factors not included in f_total

yield_ = 0.8
f_weight = 1.0

# Efficiency factors already applied in the simulation must
# not be double-counted in f_total

f_season = {
    "season" : 0.75,
    "break" : 0.25,  # Using 0.05 would double-count the weather cuts
}
f_field = {}

# Weather cut is approximately included in the sim but we want to remove
# it and replace with the factor included in f_total

f_weather_sim = {}

# LAT wide f_total from
# https://docs.google.com/spreadsheets/d/1QtA3dp0GGRWWO1WEGKVnvEYFx-MhwpxRN3enGbCDOi8/edit?pli=1&gid=845204240#gid=845204240
# on 12/10/2024

f_sensitivity = yield_ * f_weight
f_total["lat_wide"] = {
    "season" : {
        "f020" : 0.250 / f_sensitivity,
        "f030" : 0.250 / f_sensitivity,
        "f040" : 0.250 / f_sensitivity,
        "f090" : 0.250 / f_sensitivity,
        "f150" : 0.250 / f_sensitivity,
        "f220" : 0.218 / f_sensitivity,
        "f280" : 0.150 / f_sensitivity,
    },
    "break" : {
        "f020" : 0.022 / f_sensitivity,
        "f030" : 0.022 / f_sensitivity,
        "f040" : 0.022 / f_sensitivity,
        "f090" : 0.022 / f_sensitivity,
        "f150" : 0.022 / f_sensitivity,
        "f220" : 0.017 / f_sensitivity,
        "f280" : 0.012 / f_sensitivity,
    },
}
f_field["lat_wide"] = {
    "season" : 0.975,
    "break" : 0.999,
}
fskies["lat_wide"] = 0.6

# LAT delensing_max f_total from
# https://docs.google.com/spreadsheets/d/1n2NyRSKN9OZRtLJp6FTTG66upSJUAcfMDarWwU9IYb0/edit?pli=1&gid=516713372#gid=516713372
# on 2025/01/20

f_sensitivity = yield_ * f_weight
f_total["lat_delensing_max"] = {
    "season" : {
        "f020" : 0.198 / f_sensitivity,
        "f030" : 0.198 / f_sensitivity,
        "f040" : 0.198 / f_sensitivity,
        "f090" : 0.198 / f_sensitivity,
        "f150" : 0.198 / f_sensitivity,
        "f220" : 0.172 / f_sensitivity,
        "f280" : 0.119 / f_sensitivity,
    },
    "break" : {
        "f020" : 0.019 / f_sensitivity,
        "f030" : 0.019 / f_sensitivity,
        "f040" : 0.019 / f_sensitivity,
        "f090" : 0.019 / f_sensitivity,
        "f150" : 0.019 / f_sensitivity,
        "f220" : 0.014 / f_sensitivity,
        "f280" : 0.009 / f_sensitivity,
    },
}
f_field["lat_delensing_max"] = {
    "season" : 0.803,
    "break" : 0.835,
}
fskies["lat_delensing_max"] = 0.03

# f_weather from prune_schedule.py

cut_3mm_season = 5418 / 6642
cut_2mm_season = 4568 / 6642
cut_3mm_break = 474 / 2234
cut_2mm_break = 401 / 2234
f_weather_sim["lat_wide"] = {
    "season" : {
        "f020" : cut_3mm_season,
        "f030" : cut_3mm_season,
        "f040" : cut_3mm_season,
        "f090" : cut_3mm_season,
        "f150" : cut_3mm_season,
        "f220" : cut_2mm_season,
        "f280" : cut_2mm_season,
    },
    "break" : {
        "f020" : cut_3mm_break,
        "f030" : cut_3mm_break,
        "f040" : cut_3mm_break,
        "f090" : cut_3mm_break,
        "f150" : cut_3mm_break,
        "f220" : cut_2mm_break,
        "f280" : cut_2mm_break,
    },
}

cut_3mm_season = 13011 / 15909
cut_2mm_season = 10955 / 15909
cut_3mm_break =  1158 / 5414
cut_2mm_break =  908 / 5414
f_weather_sim["lat_delensing_max"] = {
    "season" : {
        "f020" : cut_3mm_season,
        "f030" : cut_3mm_season,
        "f040" : cut_3mm_season,
        "f090" : cut_3mm_season,
        "f150" : cut_3mm_season,
        "f220" : cut_2mm_season,
        "f280" : cut_2mm_season,
    },
    "break" : {
        "f020" : cut_3mm_break,
        "f030" : cut_3mm_break,
        "f040" : cut_3mm_break,
        "f090" : cut_3mm_break,
        "f150" : cut_3mm_break,
        "f220" : cut_2mm_break,
        "f280" : cut_2mm_break,
    },
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

# Number of telescope years is a parameter

n_years = {
    "lat_wide" : [14],
    "lat_delensing_max" : [16, 26, 36],
}


# Loop over all covariance matrices

for flavor in "lat_wide", "lat_delensing_max":
    for n_year in n_years[flavor]:
        nrow, ncol = 2, 4
        fig = plt.figure(figsize=[4 * ncol, 4 * nrow])
        fig.suptitle(f"{flavor} {n_year} years")
        iplot = 0
        for band in f_total[flavor]["season"].keys():
            covs = []
            for period in "season", "break":
                fname_in = f"outputs/{flavor}/{band}/{period}/mapmaker_cov.fits"
                print(f"Loading {fname_in}")
                cov = hp.read_map(fname_in, [0, 3, 5])  # II, QQ, UU

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
        fname_plot = f"plots/{flavor}_{n_year}years.png"
        fig.savefig(fname_plot)
        print(f"Plot saved in {fname_plot}")
