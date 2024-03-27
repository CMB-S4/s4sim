import glob
import os
import pickle
import sys

from toast.ops.polyfilter.kernels import filter_polynomial


fnames = sorted(glob.glob("filter_data*pck"))
nfile = len(fnames)

for ifile, fname in enumerate(fnames):
    print(f"Loading {ifile + 1} / {nfile} : {fname}", flush=True)
    with open(fname, "rb") as f:
        [
            order, flags, signals, starts, stops
        ] = pickle.load(f)
    print(f"Filtering {fname} ... ", end="", flush=True)
    filter_polynomial(order, flags, signals, starts, stops)
    print(f" done!", flush=True)
