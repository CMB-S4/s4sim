import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np

hs = hp.read_map("../pole_sat/out/00000000/toast_telescope_all_time_all_hmap.fits")
hl = hp.read_map("out/00000000/toast_telescope_all_time_all_hmap.fits")

sgood = hs > 0
nsgood = np.sum(sgood)
ordered = np.sort(hs[sgood])
lim = ordered[nsgood // 10]
mask = hs > lim
nmask = np.sum(mask)

lgood = hl > 0
nlgood = np.sum(lgood)

print()
print("Raw fsky, SAT : {}".format(nsgood / hs.size))
print("Raw fsky, LAT : {}".format(nlgood / hl.size))
print("Fraction of raw SAT not seen by LAT : {}".format(
    1 - np.sum(np.logical_and(sgood, lgood)) / nsgood
))
print("Fraction of raw LAT not seen by SAT : {}".format(
    1 - np.sum(np.logical_and(sgood, lgood)) / nlgood
))

print("Best 90% fsky, SAT : {}".format(nmask / hs.size))
print("Fraction of 90% SAT not seen by LAT : {}".format(
    1 - np.sum(np.logical_and(mask, lgood)) / nmask
))
print()

nrow, ncol = 2, 3
rot = [40, -55]
xsize = 800
reso = 8

fig = plt.figure(figsize=[8 * ncol, 6 * nrow])

m = np.copy(hs)
m[m==0] = hp.UNSEEN
hp.gnomview(m, rot=rot, xsize=xsize, reso=reso, title="SAT hits", sub=[nrow, ncol, 1])

hp.gnomview(mask, rot=rot, xsize=xsize, reso=reso, title="SAT 90% mask", sub=[nrow, ncol, 2])

m = np.logical_and(hs != 0, hl == 0)
hp.gnomview(m, rot=rot, xsize=xsize, reso=reso, title="Pixels seen only by SAT", sub=[nrow, ncol, 3])

m = np.copy(hl)
m[m==0] = hp.UNSEEN
hp.gnomview(m, rot=rot, xsize=xsize, reso=reso, title="LAT hits", sub=[nrow, ncol, 4])

m = np.logical_and(hs == 0, hl != 0)
hp.gnomview(m, rot=rot, xsize=xsize, reso=reso, title="Pixels seen only by LAT", sub=[nrow, ncol, 5])

m = np.logical_and(mask != 0, hl == 0)
hp.gnomview(m, rot=rot, xsize=xsize, reso=reso, title="90% SAT pixels not seen by LAT", sub=[nrow, ncol, 6])

fig.savefig("LAT_vs_SAT_hits.png")
