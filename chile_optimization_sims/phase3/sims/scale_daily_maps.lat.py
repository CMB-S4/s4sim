# This script scales the simulated noise covariances with efficiency and
# survey factors that were not included in the simulation

import glob
import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np

from toast.pixels_io_healpix import read_healpix


# Total target efficiency factor

f_total = {}

# fsky for measurement requirement

# fskies = {}

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
        "f090" : 0.250 / f_sensitivity,
        "f150" : 0.250 / f_sensitivity,
    },
    "break" : {
        "f090" : 0.022 / f_sensitivity,
        "f150" : 0.022 / f_sensitivity,
    },
}
f_field["lat_wide"] = {
    "season" : 0.955,
    "break" : 0.996,
}
# fskies["lat_wide"] = 0.7

# LAT delensing_max f_total from
# https://docs.google.com/spreadsheets/d/1n2NyRSKN9OZRtLJp6FTTG66upSJUAcfMDarWwU9IYb0/edit?pli=1&gid=516713372#gid=516713372
# on 2025/01/20

f_sensitivity = yield_ * f_weight
f_total["lat_delensing_max"] = {
    "season" : {
        "f090" : 0.198 / f_sensitivity,
        "f150" : 0.198 / f_sensitivity,
    },
    "break" : {
        "f090" : 0.019 / f_sensitivity,
        "f150" : 0.019 / f_sensitivity,
    },
}
f_field["lat_delensing_max"] = {
    "season" : 0.803,
    "break" : 0.835,
}

# No specific calculations for supplemental fields

f_total["lat_delensing_sun90bk"] = f_total["lat_delensing_max"]
f_total["lat_wide_supplement"] = f_total["lat_wide"]
f_total["lat_roman_supplement"] = f_total["lat_delensing_max"]

f_field["lat_wide_supplement"] = f_field["lat_wide"]

f_field["lat_delensing_sun90bk"] = {
    "season" : 0.624,
    "break" : 0.623,
}

f_field["lat_roman_supplement"] = {
    "season" : 0.426,
    "break" : 0.523,
}

# fskies["lat_delensing_max"] = 0.03

# f_weather from prune_schedule.py

cut_3mm_season = 5371 / 6588
cut_2mm_season = 4505 / 6588
cut_3mm_break = 471 / 2230
cut_2mm_break = 400 / 2230
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

f_weather_sim["lat_wide_supplement"] = f_weather_sim["lat_wide"]

cut_3mm_season = 10096 / 12356
cut_2mm_season = 8468 / 12356
cut_3mm_break =  863 / 4036
cut_2mm_break =  634 / 4036
f_weather_sim["lat_delensing_sun90bk"] = {
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

cut_3mm_season = 3245 / 4040
cut_2mm_season = 2716 / 4040
cut_3mm_break =  407 / 1646
cut_2mm_break =  304 / 1646
f_weather_sim["lat_roman_supplement"] = {
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
    "lat_wide" : [10],
    "lat_delensing_sun90bk" : [10],
    "lat_wide_supplement" : [10],
    "lat_roman_supplement" : [10],
}


# Loop over all covariance matrices

master =  "lat_delensing_sun90bk"
#for flavor in "lat_wide", "lat_delensing_max":
for supplement in "lat_wide_supplement", "lat_roman_supplement":
    for day in sorted(glob.glob("daily_outputs/lat_wide/f090/*")):
        day = os.path.basename(day)
        if day.split("-")[1] in ["01", "02", "03"]:
            period = "break"
        else:
            period = "season"
        for band in f_total[master]["season"].keys():
            covs = []
            for flavor in master, supplement:
                fname_in = f"daily_outputs/{flavor}/{band}/{day}/mapmaker_cov.h5"
                if not os.path.isfile(fname_in):
                    print(f"No covariance: {fname_in}")
                    continue
                print(f"Loading {fname_in}")
                cov = read_healpix(fname_in)  # II

                # Scale the white noise covariance

                scale = 1.
                # Compensate for focalplane decimation
                scale /= thinfp[band]
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

            if len(covs) == 0:
                continue

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

            outdir = f"scaled_daily_outputs/delensing_and_{supplement}/{band}"
            os.makedirs(outdir, exist_ok=True)
            """
            fname_out = f"{outdir}/{day}_cov.fits"
            print(f"Writing {fname_out}")
            hp.write_map(fname_out, cov, dtype=np.float32, coord="C", overwrite=True)
            """

            # Derive depth from the covariance and save

            nside = hp.get_nside(cov)
            pix_area = hp.nside2pixarea(nside, degrees=True) * 3600  # arcmin^2
            depth = np.sqrt(cov * pix_area) * 1e6  # uK.arcmin
            fname_out = f"{outdir}/{day}_depth.fits"
            print(f"Writing {fname_out}")
            hp.write_map(
                fname_out,
                depth,
                dtype=np.float32,
                coord="C",
                overwrite=True,
            )

            plotdir = f"frames/delensing_and_{supplement}"
            os.makedirs(plotdir, exist_ok=True)
            fname_plot = f"{plotdir}/frame_{day}.png"
            depth[depth == 0] = hp.UNSEEN
            hp.mollview(
                depth,
                min=0,
                max=100,
                title=f"{day}",
                cmap="inferno",
                unit="$\mu$K.arcmin",
                xsize=1600,
                # format="%.3f",
            )
            plt.savefig(fname_plot)
            print(f"Wrote {fname_plot}")
            plt.close()
