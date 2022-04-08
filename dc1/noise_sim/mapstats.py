import os
import sys

import numpy as np
from toast.pixels_io import read_healpix
from toast.mpi import MPI

comm = MPI.COMM_WORLD
rank = comm.rank
ntask = comm.size

entries = []

fname_out = "mapstats.txt"
fnames = sys.argv[1:]
#fnames = ["/global/cfs/cdirs/cmbs4/dc/dc1/staging/noise_sim/outputs/LAT0_CHLAT/f090/RISING_SCAN_40-105-4/mapmaker_RISING_SCAN_40-105-4_noiseweighted_map.h5"]

my_fnames = fnames[rank::ntask]
for i, fname in enumerate(my_fnames):
    obs = fname.split("/")[-2]

    print(f"{rank:04} : Reading {i + 1:04} / {len(my_fnames):04} : {fname}", flush=True)
    m = read_healpix(fname, nest=True)
    w = read_healpix(fname.replace("noiseweighted_map", "invcov"), nest=True)
    good = m[0] != 0
    rms_i = np.std(m[0][good] / w[0][good])
    rms_qu = np.sqrt(np.var(m[1][good] / w[3][good]) + np.var(m[2][good] / w[5][good]))

    entries.append(f"{obs} {rms_i * 1e6} {rms_qu * 1e6} {fname}")

if rank == 0:
    print("Gathering entries", flush=True)

all_entries = comm.gather(entries)
if rank == 0:
    with open(fname_out, "w") as stats:
        for entries in all_entries:
            for entry in entries:
                stats.write(entry + "\n")
    print(f"Wrote {fname_out}", flush=True)
