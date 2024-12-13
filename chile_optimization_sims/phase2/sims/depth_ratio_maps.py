import glob
import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np

indir1 = "/global/cfs/cdirs/cmbs4/chile_optimization/simulations/phase1/noise_depth"
indir2 = "/global/cfs/cdirs/cmbs4/chile_optimization/simulations/phase2/noise_depth"
outdir = "depth_ratios"
os.makedirs(outdir, exist_ok=True)

# SAT

bands = ["f030", "f040", "f085", "f095", "f145", "f155", "f220", "f280"]
for band in bands:
    fname1 = glob.glob(f"{indir1}/sat_{band}_*years_depth.fits")[0]
    fname2 = glob.glob(f"{indir2}/sun90_{band}_*years_depth.fits")[0]
    m1 = hp.read_map(fname1)
    m2 = hp.read_map(fname2)
    ratio = np.zeros_like(m1)
    good = np.logical_and(m1 != 0, m2 != 0)
    ratio[good] = m2[good] / m1[good]
    outfile = os.path.join(
        outdir, os.path.basename(fname1).replace(".fits", "_ratio.fits")
    )
    hp.write_map(outfile, ratio, dtype=np.float32, overwrite=True)
    print(f"Wrote {outfile}")

# LAT - need to combine the two surveys
bands = ["f020", "f030", "f040", "f090", "f150", "f220", "f280"]
def depth_to_invcov(depth):
    invcov = np.zeros_like(depth)
    good = depth != 0
    invcov[good] = 1 / depth[good]**2
    return invcov
def invcov_to_depth(invcov):
    depth = np.zeros_like(invcov)
    good = invcov != 0
    depth[good] = 1 / invcov[good]**.5
    return depth

for band in bands:
    fname1 = f"{indir1}/lat_wide_{band}_14years_depth.fits"
    fname2 = f"{indir2}/lat_wide_{band}_14years_depth.fits"
    wide1 = hp.read_map(fname1)
    wide2 = hp.read_map(fname2)
    iwide1 = depth_to_invcov(wide1)
    iwide2 = depth_to_invcov(wide2)
    for nyear in 16, 26, 36:
        fname1 = f"{indir1}/lat_delensing_{band}_{nyear}years_depth.fits"
        fname2 = f"{indir2}/lat_delensing_{band}_{nyear}years_depth.fits"
        delens1 = hp.read_map(fname1)
        delens2 = hp.read_map(fname2)
        idelens1 = depth_to_invcov(delens1)
        idelens2 = depth_to_invcov(delens2)
        m1 = iwide1 + idelens1
        m2 = iwide2 + idelens2
        m1 = invcov_to_depth(m1)
        m2 = invcov_to_depth(m2)
        ratio = np.zeros_like(m1)
        good = np.logical_and(m1 != 0, m2 != 0)
        ratio[good] = m2[good] / m1[good]
        outfile = os.path.join(
            outdir, os.path.basename(fname1).replace(".fits", "_ratio.fits")
        )
        hp.write_map(outfile, ratio, dtype=np.float32, overwrite=True)
        print(f"Wrote {outfile}")
