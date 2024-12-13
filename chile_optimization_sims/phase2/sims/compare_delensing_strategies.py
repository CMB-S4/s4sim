import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np


h0 = hp.read_map("outputs/sun90/f145/season/mapmaker_hits.fits") * 1.
nside = hp.get_nside(h0)
pix0 = np.argmax(h0)
vec0 = hp.pix2vec(nside, pix0)
pix_core = hp.query_disc(nside, vec0, np.radians(5))
pix_wide = hp.query_disc(nside, vec0, np.radians(20))
mask_core = np.zeros(h0.size, dtype=bool)
mask_wide = np.zeros(h0.size, dtype=bool)
mask_core[pix_core] = True
mask_wide[pix_wide] = True
h0 /= np.amax(h0)
h0[h0 == 0] = hp.UNSEEN

h1 = hp.read_map("outputs/lat_delensing/f150/season/mapmaker_hits.fits") * 1.
h1 *= 0.95  # Account for turnarounds
nhit1_core = np.sum(h1[mask_core])
nhit1_wide = np.sum(h1[mask_wide])
h1norm = np.amax(h1)
h1 /= h1norm
h1[h1 == 0] = hp.UNSEEN

h2 = hp.read_map("outputs/lat_delensing_core/f150/season/mapmaker_hits.fits") * 1.
h2 *= 0.93
nhit2_core = np.sum(h2[mask_core])
nhit2_wide = np.sum(h2[mask_wide])
h2norm = np.amax(h2)
h2 /= h2norm
h2[h2 == 0] = hp.UNSEEN

h3 = hp.read_map("outputs/lat_delensing_tiled/f150/season/mapmaker_hits.fits") * 1.
h3 *= 0.91
nhit3_core = np.sum(h3[mask_core])
nhit3_wide = np.sum(h3[mask_wide])
h3norm = np.amax(h3)
h3 /= h3norm
h3[h3 == 0] = hp.UNSEEN

nrow, ncol = 2, 4
fig = plt.figure(figsize=[6 * ncol, 4 * nrow])
hp.mollview(h0, min=0, max=1, sub=[nrow, ncol, 1], title="SAT 145GHz relative hits")
hp.mollview(h1, min=0, max=1, sub=[nrow, ncol, 2], title=f"delens sums : {nhit1_core:.2e}, {nhit1_wide:.2e}")
hp.mollview(h2, min=0, max=1, sub=[nrow, ncol, 3], title=f"delens_core sums : {nhit2_core:.2e}, {nhit2_wide:.2e}")
hp.mollview(h3, min=0, max=1, sub=[nrow, ncol, 4], title=f"delens_tiled sums : {nhit3_core:.2e}, {nhit3_wide:.2e}")

def mapdiff(m1, m2):
    m1 = m1.copy()
    m1[m1 == hp.UNSEEN] = 0
    m2 = m2.copy()
    m2[m2 == hp.UNSEEN] = 0
    dmap = m1 - m2
    dmap[np.logical_and(m1 == 0, m2 == 0)] = hp.UNSEEN
    return dmap

hp.mollview(h0, min=0, max=1, sub=[nrow, ncol, 5])
hp.mollview(mask_core*1. + mask_wide*1., sub=[nrow, ncol, 5], title="Masks", reuse_axes=True, alpha=np.ones_like(h0)*.25, cbar=False)
amp = 0.3
hp.mollview(mapdiff(h1, h0), min=-amp, max=amp, sub=[nrow, ncol, 6], title="delens - SAT")
hp.mollview(mapdiff(h2, h0), min=-amp, max=amp, sub=[nrow, ncol, 7], title="delens_core - SAT")
hp.mollview(mapdiff(h3, h0), min=-amp, max=amp, sub=[nrow, ncol, 8], title="delens_tiled - SAT")
