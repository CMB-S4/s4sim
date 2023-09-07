import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np

from toast.pixels_io_healpix import read_healpix


#fname_schedule = "../scan_strategy/chile_lat/schedules/chile_schedule_lat.pruned.txt"
fname_schedule = "../scan_strategy/chile_lat/schedules/chile_schedule_lat.txt"


# /global/cfs/cdirs/cmbs4/dc/dc1/staging/noise_sim/outputs_rk/LAT0_CHLAT/f150/RISING_SCAN_40-231-12/mapmaker_RISING_SCAN_40-231-12_hits.h5
indir = "/global/cfs/cdirs/cmbs4/dc/dc1/staging/noise_sim/outputs_rk/LAT0_CHLAT/f150"
#indir = "/global/cfs/cdirs/cmbs4/dc/dc1/staging/foreground_sim/outputs_rk/LAT0_CHLAT/f150"

outdir = "hitmaps"
os.mkdirs(outdir, exist_ok=True)


schedule = []
with open(fname_schedule, "r") as file_in:
    for iline, line in enumerate(file_in):
        if iline < 3:
            continue
        parts = line.split()
        start_date = parts[0]
        start_time = parts[1]
        start = f"{start_date} {start_time}"
        target = parts[7]
        scan = parts[21]
        subscan = parts[22]
        observation = f"{target}-{scan}-{subscan}"
        schedule.append((start, observation))

n = 24
nside = 4096
npix = 12 * nside ** 2
hitmaps = []
total = np.zeros(npix, dtype=np.int32)

fig = plt.figure(figsize=[12, 8])

for start, observation in schedule[3500:4100]:
    fname_hitmap = os.path.join(
        indir,
        observation,
        f"mapmaker_{observation}_hits.h5"
    )
    if len(hitmaps) == n:
        oldest = hitmaps.pop(0)
        if oldest is not None:
            total -= oldest
    if os.path.isfile(fname_hitmap):
        print(f"Reading {fname_hitmap}")
        hitmap = read_healpix(fname_hitmap, dtype=np.int32, nest=True)[0]
        total += hitmap
    else:
        print(f"No hits for {start}")
        hitmap = None
    hitmaps.append(hitmap)
    fig.clf()
    hp.mollview(
        total,
        fig=fig.number,
        title=start,
        nest=True,
        cmap="magma",
        max=15000,
        xsize=1600,
    )
    fname_plot = f"{outdir}/hits_{start}.png"
    fig.savefig(fname_plot)
    print(f"Wrote {fname_plot}")
