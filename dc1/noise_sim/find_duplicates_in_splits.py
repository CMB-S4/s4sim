import os
import sys

from mpi4py import MPI
import numpy as np

from toast.pixels_io_healpix import read_healpix


comm = MPI.COMM_WORLD
rank = comm.rank
ntask = comm.size


#fileroot = "/global/cfs/cdirs/cmbs4/dc/dc1/staging/noise_sim/"
fileroot = None

path1 = set()
basename1 = set()
with open("coadd_maps_f150_1of2.txt", "r") as f:
    for fname in f.readlines():
        fname = fname.strip()
        if fname in path1:
            raise RuntimeError(f"{fname} is already listed!")
        path1.add(fname)
        fname = os.path.basename(fname)
        if fname in basename1:
            raise RuntimeError(f"{fname} is already listed!")
        basename1.add(fname)

path2 = set()
basename2 = set()
with open("coadd_maps_f150_2of2.txt", "r") as f:
    for fname in f.readlines():
        fname = fname.strip()
        if fname in path1:
            raise RuntimeError(f"{fname} is already listed!")
        path2.add(fname)
        fname = os.path.basename(fname)
        if fname in basename1:
            raise RuntimeError(f"{fname} is already listed!")
        basename2.add(fname)

f_overlap = open(f"with_overlap.{rank}.txt", "w")
f_no_overlap = open(f"wo_overlap.{rank}.txt", "w")

i = -1
for fname1 in sorted(path1):
    print(f"{rank:4} : fname1 = {fname1}")
    m1 = read_healpix(fname1, field=[0], nest=True, verbose=False)
    if fileroot is not None:
        fname1 = fname1.replace(fileroot, "")
    good1 = m1 != 0
    ngood1 = np.sum(good1)
    for fname2 in sorted(path2):
        i += 1
        if i % ntask != rank:
            continue
        print(f"{rank:4} :   fname2 = {fname2}")
        m2 = read_healpix(fname2, field=[0], nest=True, verbose=False)
        if fileroot is not None:
            fname2 = fname2.replace(fileroot, "")
        good2 = m2 != 0
        ngood2 = np.sum(good2)
        good12 = good1 * good2
        ngood12 = np.sum(good12)
        frac = ngood12 / max(ngood1, ngood2)
        if frac > 0.25:
            rms_sum = np.std((m1 + m2)[good12])
            rms_diff = np.std((m1 - m2)[good12])
            f_overlap.write(
                f"{fname1} {fname2} {frac} {rms_sum} {rms_diff} "
                f"{rms_diff / rms_sum}\n"
            )
            f_overlap.flush()
        else:
            f_no_overlap.write(f"{fname1} {fname2}\n")
            f_no_overlap.flush()
f_overlap.close()
f_no_overlap.close()
