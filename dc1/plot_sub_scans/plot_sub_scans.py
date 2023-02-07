import os
import sys

import astropy.units as u
import healpy as hp
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
times -= times[0]

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
os.makedirs("frames", exist_ok=True)
while tstart < times[-1]:
    tstop = tstart + time_step
    istart, istop = np.searchsorted(times, [tstart, tstop])
    ind = slice(istart, istop)
    temp_data = trim_data(data, ind)
    pixel_pointing.exec(temp_data)
    for obs in temp_data.obs:
        pixels = obs.detdata["pixels"].data.ravel()
        good = pixels >= 0
        bin_hits(hitmap, pixels[good].copy())

    fname = f"frames/frame_{frame:04}.png"
    plt.clf()
    hp.mollview(hitmap, title=f"t = {tstart:5.0f} s", min=0, max=hit_max, cmap="magma")
    plt.savefig(fname)
    print(f"Wrote {fname}, {tstart / times[-1] * 100.:6.1f}% done", flush=True)
    
    tstart = tstop
    frame += 1

print("Done. Now animate with\n  convert -loop 0 -delay 10 frames/frame*png movie.gif")
