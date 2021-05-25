#!/usr/bin/env python

from glob import glob
import os
import sys
from time import time

import numpy as np
import scipy.sparse


if len(sys.argv) != 2:
    print("Usage: combine_observation_matrix.py <rootname>")
    sys.exit()

t0 = time()

rootname = sys.argv[1]

datafiles = sorted(glob(f"{rootname}.*.*.*.data.npy"))

all_data = []
all_indices = []
all_indptr = [0]

current_row = 0
current_offset = 0
shape = None

t1 = time()

for datafile in datafiles:
    parts = datafile.split(".")
    row_start = int(parts[-5])
    row_stop = int(parts[-4])
    nrow_tot = int(parts[-3])
    if shape is None:
        shape = (nrow_tot, nrow_tot)
    elif shape[0] != nrow_tot:
        raise RuntimeError("Mismatch in shape")
    if current_row != row_start:
        all_indptr.append(np.zeros(row_start - current_row) + current_offset)
        current_row = row_start
    print(f"Loading {datafile}", flush=True)
    data = np.load(datafile)
    indices = np.load(datafile.replace(".data.", ".indices."))
    indptr = np.load(datafile.replace(".data.", ".indptr."))
    all_data.append(data)
    all_indices.append(indices)
    indptr += current_offset
    all_indptr.append(indptr[1:])
    current_row = row_stop
    current_offset = indptr[-1]

print(f"Parts loaded in {time() - t1:.1f} s", flush=True)

if current_row != nrow_tot:
    all_indptr.append(np.zeros(nrow_tot - current_row) + current_offset)

t1 = time()
print("Constructing CSR matrix", flush=True)

all_data = np.hstack(all_data)
all_indices = np.hstack(all_indices)
all_indptr = np.hstack(all_indptr)
obs_matrix = scipy.sparse.csr_matrix((all_data, all_indices, all_indptr), shape)

print(f"Constructed in {time() - t1:.1f} s", flush=True)

t1 = time()
print(f"Writing {rootname}.npz", flush=True)
scipy.sparse.save_npz(rootname, obs_matrix)
print(f"Wrote in {time() - t1:.1f} s", flush=True)

print(f"Done in {time() - t0:.1f} s!", flush=True)
