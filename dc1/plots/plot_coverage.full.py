from glob import glob
import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from toast.pixels_io import read_healpix

overwrite = True
ort = False
tele = "chlat"
#tele = "splat"
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
print(f"Found {len(dirs)} directories matching pattern '{pattern}'", flush=True)
dirs = sorted(dirs)

if tele == "chlat":
    window = 24 * 1
    step = 1000
    last = 5700
    ncluster = 24
elif tele == "spsat":
    window = 10
    step = 1
    last = 12
    ncluster = 12
else:
    window = 10
    step = 1
    last = 100
    ncluster = 1

all_hits = None
offset = 0
obs = 0
nleft = ncluster
while True:
    indir = dirs[obs]
    print(indir, nleft, obs, flush=True)
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
        all_hits = np.zeros(hits.size, dtype=np.uint8)
    all_hits[hits != 0] = 1
    if nleft > 1:
        nleft -= 1
    else:
        obs += step - ncluster
        nleft = ncluster
    obs += 1
    if obs > last:
        break

if ort:
    fname_out = f"full_hits.{tele}.ortho.png"
else:
    fname_out = f"full_hits.{tele}.png"
if ort:
    hp.orthview(all_hits, nest=True, cmap="magma", title=None, cbar=False, rot=[0, 30, 30])
    hp.graticule(45)
else:
    hp.mollview(all_hits, nest=True, cmap="magma", title=None, cbar=False)
plt.savefig(fname_out)
print(f"Wrote {fname_out}", flush=True)
plt.close()

fname_out = f"full_hits.{tele}.fits"
hp.write_map(fname_out, all_hits, dtype=np.uint8, nest=True, coord="C", overwrite=True)
print(f"Wrote {fname_out}", flush=True)
