from glob import glob
import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from toast.pixels_io import read_healpix

from mpi4py import MPI

comm = MPI.COMM_WORLD

overwrite = True
ort = False
#tele = "chlat"
tele = "splat"
#tele = "spsat"

if tele == "spsat":
    #pattern = f"/global/cfs/cdirs/cmbs4/dc/dc1/staging/cmb_sim/outputs/{tele}/f085/{tele}_split_schedule_????"
    pattern = f"../cmb_sim/outputs/{tele}/f085/{tele}_split_schedule_????"
elif tele == "splat":
    #pattern = f"/global/cfs/cdirs/cmbs4/dc/dc1/staging/cmb_sim/outputs/{tele}/f085/{tele}_split_schedule_????"
    pattern = f"../cmb_sim/outputs/{tele}/f090/{tele}_split_schedule_????"
else:
    pattern = f"/global/cfs/cdirs/cmbs4/dc/dc1/staging/cmb_sim/outputs/{tele}/f090/{tele}_split_schedule_????"
dirs = glob(pattern)
if comm.rank == 0:
    print(f"Found {len(dirs)} directories matching pattern '{pattern}'", flush=True)
dirs = sorted(dirs)

if tele == "chlat":
    window = 24 * 1
    step = 400
elif tele == "spsat":
    window = 10
    step = 15
else:
    window = 10
    step = 10

all_hits = None
offset = comm.rank * step
obs = max(0, offset - window)
for indir in list(dirs[obs : ]):
    print(indir, flush=True)
    try:
        if tele == "spsat":
            fname_in = glob(f"{indir}/filterbin*hits.h5")[0]
        else:
            fname_in = glob(f"{indir}/mapmaker*hits.h5")[0]
    except Exception as e:
        print(f"Failed to read hits from {indir}: {e}", flush=True)
        obs += 1
        continue
    hits = read_healpix(fname_in, nest=True)[0]
    if all_hits is None:
        all_hits = np.zeros(hits.size)
    all_hits[all_hits != 0] -= 1
    all_hits[all_hits < 0] = 0
    all_hits[hits != 0] = window
    if ort:
        fname_out = f"hits_{obs:04}.{tele}.ortho.png"
    else:
        fname_out = f"hits_{obs:04}.{tele}.png"
    if obs >= offset and (not os.path.isfile(fname_out) or overwrite):
        #if tele == "splat":
        #    title = f"{obs:04}" + " {:12}".format(
        #        os.path.basename(fname_in).replace("mapmaker_", "").replace("_hits.h5", "")
        #    )
        #else:
        title = f"{obs:04}"
        if ort:
            hp.orthview(all_hits, nest=True, cmap="magma", title=title, cbar=False, rot=[0, 30, 30])
            hp.graticule(45)
        else:
            hp.mollview(all_hits, nest=True, cmap="magma", title=title, cbar=False)
        plt.savefig(fname_out)
        print(f"Wrote {fname_out}", flush=True)
        plt.close()
    obs += 1
    if obs == offset + step:
        break

comm.barrier()
