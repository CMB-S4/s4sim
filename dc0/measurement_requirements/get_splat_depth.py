from mpi4py import MPI

import os
import pickle
import sys

import astropy.units as u
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import c, h, k

from toast.pixels_io_healpix import read_healpix


comm = MPI.COMM_WORLD
ntask = comm.size
rank = comm.rank

if len(sys.argv) > 1:
    band = sys.argv[1]
else:
    #band = "f150"
    band = "f090"
fname_out = f"depths_splat_{band}.pck"
if os.path.isfile(fname_out):
    with open(fname_out, "rb") as f:
        depths = pickle.load(f)
else:
    depths = {}

TCMB = 2.72548
ifreq = int(band[1:])
freq = ifreq * u.GHz
fwhm_arcmin = {90 : 2.5, 150 : 1.6}[ifreq]
fwhm = fwhm_arcmin * u.arcmin
sigma = fwhm / np.sqrt(8 * np.log(2))
solid_angle = 2 * np.pi * sigma**2

fname_schedule = "../scan_strategy/pole_lat/schedules/pole_schedule_lat.txt"

indir = f"/global/cfs/cdirs/cmbs4/dc/dc1/staging/noise_sim/outputs_rk/LAT2_SPLAT/{band}"

schedule = {}
with open(fname_schedule, "r") as file_in:
    for iline, line in enumerate(file_in):
        if iline < 3:
            continue
        parts = line.split()
        start_date = parts[0]
        # start_time = parts[1]
        # start = f"{start_date} {start_time}"
        target = parts[7]
        if target == "CALIBRATION_BREAK":
            continue
        scan = parts[21]
        subscan = parts[22]
        observation = f"{target}-{scan}-{subscan}"
        if start_date not in schedule:
            schedule[start_date] = []
        schedule[start_date].append(observation)

nside = 4096
npix = 12 * nside ** 2
limit = int(npix * 0.03)  # Measured over best 3% of the sky
pix_area = hp.nside2pixarea(nside)
year = 2027

# Delta-function bandpass

nu_cmb = k * TCMB / h
alpha = 2 * k**3 * TCMB**2 / h**2 / c**2
x = freq.to_value(u.Hz) / nu_cmb
db_dt = alpha * x**4 * np.exp(x) / (np.exp(x) - 1) ** 2
kcmb2Jysr = db_dt * 1e26
kcmb2mJysr = kcmb2Jysr * 1e3

job = -1
for month in range(1, 13):
# for month in range(1, 2):
    for day in range(1, 32):
        job += 1
        if job % ntask != rank:
            continue
        date = f"{year}-{month:02}-{day:02}"
        print(f"{rank} : {date} processing", flush=True)
        if date in depths:
            # Already processed before
            print(f"{rank} : {date} already processed", flush=True)
            continue
        if date not in schedule:
            # not a scheduled date
            print(f"{rank} : {date} not scheduled", flush=True)
            # depths[date] = None
            continue
        daily_invcov = np.zeros(npix)
        for observation in schedule[date]:
            fname_invcov = os.path.join(
                indir,
                observation,
                f"mapmaker_{observation}_invcov.h5"
            )
            if os.path.isfile(fname_invcov):
                print(f"{rank} : Reading {fname_invcov}", flush=True)
                invcov = read_healpix(
                    fname_invcov, [0], dtype=np.float32, nest=True
                )[0]
                daily_invcov += invcov
            else:
                print(f"{rank} : Not found: {fname_invcov}", flush=True)
        good = daily_invcov != 0
        ngood = np.sum(good)
        if ngood < limit:
            depth = np.inf
        else:
            depth = np.sqrt(1 / daily_invcov[good]) * kcmb2mJysr \
                * solid_angle.to_value(u.steradian)
            # depth = np.sort(depth)[limit]
            depth = np.median(np.sort(depth)[:limit])
            print(f"{rank} : {date} : Depth = {depth:.2f} mJy", flush=True)
        depths[date] = depth

for rank_send in range(1, ntask):
    if rank_send == rank:
        comm.send(depths, dest=0, tag=rank)
    elif rank == 0:
        depths_recv = comm.recv(source=rank_send, tag=rank_send)
        depths.update(depths_recv)
    comm.barrier()

if rank == 0:
    with open(fname_out, "wb") as f:
        pickle.dump(depths, f)
    print(f"{rank} : Wrote {fname_out}", flush=True)
