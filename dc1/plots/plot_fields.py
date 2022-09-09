from glob import glob
import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np

#tele, TELE = "spsat", "SAT"
#tele, TELE = "chlat", "LAT"
tele, TELE = "splat", "LAT"

#cmap = "bwr"
#cmap = plt.get_cmap("jet").copy()
#cmap.set_bad("whitesmoke")
#cmap = hp.projaxes.create_colormap("jet", badcolor="white", bgcolor="white")
cmap = "jet"

rootdir = "/global/cfs/cdirs/cmbs4/dm/dstool_202102/input_pysm/"

if tele == "spsat":
    nside = 512
else:
    nside = 4096

field = hp.read_map(f"full_hits.{tele}.fits")

fnames = glob(f"{rootdir}/{nside}/combined_foregrounds/0000/*fits")
for fname_fg in fnames:
    print(fname_fg)
    basename = os.path.basename(fname_fg)
    band = basename.split("_")[4].split("-")[1]
    if tele == "splat" and "P" not in band:
        continue
    if tele == "chlat" and "P" in band:
        continue
    print(band)
    m = hp.read_map(fname_fg, None)
    freq = {
        "ULFL1" : 20,
        "LFL1" : 25,
        "LFL2" : 40,
        "MFL1" : 90,
        "MFL2" : 150,
        "HFL1" : 230,
        "HFL2" : 280,
    }[band.replace("P", "")]
    print(freq)
    fname_cmb = f"{rootdir}/{nside}/cmb/0000/cmbs4_cmb_uKCMB_{TELE}-{band}_nside{nside}_0000.fits"
    #fname_cmb = f"{rootdir}/{nside}/cmb_unlensed_solardipole/0000/cmbs4_cmb_unlensed_solardipole_uKCMB_SAT-{band}_nside{nside}_0000.fits"
    print(fname_cmb)
    m += hp.read_map(fname_cmb, None)
    amp = 300
    #frac = 0.0
    #alpha = (frac + field) / (frac + 1)
    alpha = None
    m[:, field == 0] = hp.UNSEEN

    common = {
        "alpha" : alpha,
        "title" : None,
        "cmap" : cmap,
        "cbar" : False,
        "hold" : True,
        "norm" : "hist",
        "xsize" : 1600,
        "badcolor" : "gainsboro",
        "bgcolor" : "white",
    }

    figsize = [12, 8]
    fig = plt.figure(figsize=figsize)
    hp.mollview(m[0], **common)
    fig.text(0.39, 0.65, f"{freq:03}GHz I", size=42)
    fig.savefig(f"field_{tele}_{freq:03}.I.png", transparent=True)
    plt.close()

    amp = 10
    fig = plt.figure(figsize=figsize)
    hp.mollview(m[1], **common)
    fig.text(0.39, 0.65, f"{freq:03}GHz Q", size=42)
    fig.savefig(f"field_{tele}_{freq:03}.Q.png", transparent=True)
    plt.close()

    fig = plt.figure(figsize=figsize)
    hp.mollview(m[2], **common)
    fig.text(0.39, 0.65, f"{freq:03}GHz U", size=42)
    fig.savefig(f"field_{tele}_{freq:03}.U.png", transparent=True)
    plt.close()

# Make an empty plot for CHLAT 20GHz

common = {
    "alpha" : alpha,
    "title" : None,
    "cmap" : cmap,
    "cbar" : False,
    "hold" : True,
    "norm" : "hist",
    "xsize" : 1600,
    "badcolor" : "gainsboro",
    "bgcolor" : "white",
}

tele = "chlat"
freq = 20

figsize = [12, 8]
fig = plt.figure(figsize=figsize)
hp.mollview(m[0], **common)
fig.text(0.39, 0.65, f"{freq:03}GHz I", size=42)
fig.savefig(f"field_{tele}_{freq:03}.I.png", transparent=True)
plt.close()

amp = 10
fig = plt.figure(figsize=figsize)
hp.mollview(m[1], **common)
fig.text(0.39, 0.65, f"{freq:03}GHz Q", size=42)
fig.savefig(f"field_{tele}_{freq:03}.Q.png", transparent=True)
plt.close()

fig = plt.figure(figsize=figsize)
hp.mollview(m[2], **common)
fig.text(0.39, 0.65, f"{freq:03}GHz U", size=42)
fig.savefig(f"field_{tele}_{freq:03}.U.png", transparent=True)
plt.close()

