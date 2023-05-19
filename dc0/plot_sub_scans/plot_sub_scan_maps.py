import os
import sys

import astropy.units as u
import healpy as hp
import matplotlib
import matplotlib.pyplot as plt
from mpi4py import MPI
import numpy as np
# Numba is not available for Python 3.11 ...
from numba import jit

import toast
from toast.ops import LoadHDF5


comm = MPI.COMM_WORLD
rank = comm.rank
ntask = comm.size

inroot = "/global/cfs/cdirs/cmbs4/dc/dc1/staging/noise_sim/outputs_rk/LAT0_CHLAT/f150/"
indir = "SETTING_SCAN_40-156-12"
inpath = os.path.join(inroot, indir)
data = toast.Data(toast.mpi.Comm(world=MPI.COMM_SELF))

nside = 512
nest = False

matplotlib.rcParams["font.family"] = "monospace"

print(f"Loading data from {inpath}", flush=True)
loader = LoadHDF5(
    volume=inpath,
    detdata=["noexistent entry"],
    shared=["times", "flags", "boresight_radec"],
    pattern="^(RISING|SETTING).*h5",
)
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
        obs_out.shared.create_column(
            "boresight_radec", shape=(nsample, 4), dtype=quats.dtype
        )
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

@jit(nopython=True)
def bin_hits(hitmap, pixels):
    for pixel in pixels:
        hitmap[pixel] += 1

npix = 12 * nside**2
hitmap = np.zeros(npix)

print(f"Plotting hits", flush=True)

tstart = times[0]
framedir = f"sub_scans_{indir}"
os.makedirs(framedir, exist_ok=True)
intervals = data.obs[0].intervals["throw"]
ninterval = len(intervals)
for frame, interval in enumerate(intervals):
    if frame % ntask != rank:
        continue
    istart = interval.first
    istop = interval.last + 1
    tstart = interval.start - time_start
    tstop = interval.stop - time_start
    ind = slice(istart, istop)
    temp_data = trim_data(data, ind)
    pixel_pointing.exec(temp_data)
    hitmap[:] = 0
    for obs in temp_data.obs:
        pixels = obs.detdata["pixels"].data.ravel()
        good = pixels >= 0
        bin_hits(hitmap, pixels[good].copy())
    del temp_data

    fname = f"{framedir}/sub_scan_{frame:04}.png"
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
        cmap="magma",
        unit="hits",
        lonra=[-30, 30],
        latra=[-70, 30],
    )
    plt.savefig(fname)
    plt.close()
    print(f"Wrote {fname}, {tstart / times[-1] * 100.:6.1f}% done", flush=True)

comm.Barrier()
if rank == 0:
    print("Done.")
