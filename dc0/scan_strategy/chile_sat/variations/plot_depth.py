import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np


m = hp.read_map("outputs/baseline/mapmaker_cov.fits")
good = m != 0
m[good] = 1 / m[good]
m[good] /= np.amax(m[good])
m[m==0] = hp.UNSEEN

bk = hp.read_map("/home/reijo/work/bicep/bk18_mask_largefield_cel_n0512.fits")

hp.mollview(
    m,
    cbar=False,
    notext=True,
    title="",
    cmap="inferno",
)
plt.savefig("baseline_depth.png")

hp.mollview(
    m,
    alpha=np.ones_like(m)*.5,
    cbar=False,
    notext=True,
    title="",
    cmap="inferno",
)

hp.mollview(
    bk,
    reuse_axes=True,
    cbar=False,
    alpha=np.isfinite(bk) * .8,
    notext=True,
    title="",
    cmap="inferno",
)
plt.savefig("baseline_depth+bicep.png")
