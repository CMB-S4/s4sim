import os
import sys

import astropy.units as u
import healpy as hp
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
# Numba is not available for Python 3.11 ...
# from numba import jit

import toast
from toast.ops import LoadHDF5

indir = "./data"
data = toast.Data()

nside = 512
nest = False
time_step = 5.0

matplotlib.rcParams["font.family"] = "monospace"

print(f"Loading data from {indir}", flush=True)
loader = LoadHDF5(volume=indir, detdata=[], shared=["times", "flags", "boresight_radec"])
loader.exec(data)

def trim_data(data_in, ind):
    data_out = toast.Data()
    for obs_in in data_in.obs:
        times = obs_in.shared["times"].data[ind]
        flags = obs_in.shared["flags"].data[ind]
        quats = obs_in.shared["boresight_radec"].data[ind]
        nsample = times.size
        obs_out = toast.Observation(
            obs_in.comm,
            obs_in.telescope,
            nsample,
            name=obs_in.name,
            uid=obs_in.uid,
        )
        obs_out.shared.create_column("times", shape=(nsample,), dtype=times.dtype)
        obs_out.shared["times"].set(times, offset=(0,), fromrank=0)
        obs_out.shared.create_column("flags", shape=(nsample,), dtype=flags.dtype)
        obs_out.shared["flags"].set(flags, offset=(0,), fromrank=0)
        obs_out.shared.create_column("boresight_radec", shape=(nsample, 4), dtype=quats.dtype)
        obs_out.shared["boresight_radec"].set(quats, offset=(0, 0), fromrank=0)
        data_out.obs.append(obs_out)
    return data_out

quat_pointing = toast.ops.PointingDetectorSimple()
pixel_pointing = toast.ops.PixelsHealpix(
    detector_pointing=quat_pointing, nside=nside, nest=nest
)

obs = data.obs[0]
times = obs.shared["times"]
time_start = times[0]
times -= time_start

# @jit(nopython=True)
def bin_hits(hitmap, pixels):
    for pixel in pixels:
        hitmap[pixel] += 1

npix = 12 * nside**2
hitmap = np.zeros(npix)

print(f"Plotting hits", flush=True)

#hitmap_tot = np.zeros(npix)
#pixels = obs.detdata["pixels"].data.ravel()
#good = pixels >= 0
#bin_hits(hitmap_tot, pixels[good].copy())
#hit_max = np.amax(hitmap_tot)
hit_max = None

tstart = times[0]
frame = 0
os.makedirs("sub_scans", exist_ok=True)
intervals = data.obs[0].intervals["throw"]
ninterval = len(intervals)
for interval in intervals:
    istart = interval.first
    istop = interval.last + 1
    tstart = interval.start - time_start
    tstop = interval.stop - time_start
    ind = slice(istart, istop)
    temp_data = trim_data(data, ind)
    pixel_pointing.exec(temp_data)
    if frame % 2 == 0:
        hitmap[:] = 0
    for obs in temp_data.obs:
        pixels = obs.detdata["pixels"].data.ravel()
        good = pixels >= 0
        bin_hits(hitmap, pixels[good].copy())

    fname = f"sub_scans/sub_scan_{frame:04}.png"
    #hitmap[hitmap == 0] = hp.UNSEEN
    #hp.mollview(
    #    hitmap,
    #    title=f"sub scan {frame+1:3} / {ninterval},  t = {tstart:4.0f} - {tstop:4.0f} s",
    #    min=0,
    #    max=hit_max,
    #    cmap="magma",
    #    unit="hits",
    #)
    hp.cartview(
        hitmap,
        title=f"sub scan {frame+1:2} / {ninterval},  t = {tstart:4.0f} - {tstop:4.0f} s",
        min=0,
        max=hit_max,
        cmap="magma",
        unit="hits",
        lonra=[90, 180],
        latra=[-70, 30],
    )
    plt.savefig(fname)
    plt.close()
    print(f"Wrote {fname}, {tstart / times[-1] * 100.:6.1f}% done", flush=True)
    
    tstart = tstop
    frame += 1

print("Done.")
