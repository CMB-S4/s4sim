import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np

hits = hp.read_map("out_v1/00000000/pole_noise_SAT_MFLS1_filtered_telescope_all_time_all_hmap.fits")
m1 = hp.read_map("out-test/00000000/pole_atmosphere_SAT_MFLS1_thinfp1_telescope_all_time_all_bmap.fits",None)
m4 = hp.read_map("out-test/00000000/pole_atmosphere_SAT_MFLS1_thinfp4_telescope_all_time_all_bmap.fits",None)

amp = 100

good = hits > 0
sorted_hits = np.sort(hits[good])
ngood = sorted_hits.size
lim = sorted_hits[int(.5 * ngood)]
mask = hits > lim
rms1 = np.std(m1[0][mask]) * 1e6
rms4 = np.std(m4[0][mask]) * 1e6
rms14 = np.std((m4-m1)[0][mask]) * 1e6

nrow, ncol = 2, 2
fig = plt.figure(figsize=[8 * ncol, 6 * nrow])


args = {
    "reso" : .75,
    "xsize" : 1600,
    "rot" : [40, -50],
    "min" : -amp,
    "max" : amp,
    "unit" : "$\mu$K$_\mathrm{CMB}$",
    "cmap" : "coolwarm",
}

grargs = {
    "dpar" : 5,
    "dmer" : 15,
    "local" : False,
}

hp.gnomview(
    m1[0] * 1e6,
    **args,
    title="thinfp=1, rms={:.1f}uK".format(rms1),
    sub=[nrow, ncol, 1],
)
hp.graticule(**grargs)

hp.gnomview(
    m4[0] * 1e6,
    **args,
    title="thinfp=4, rms={:.1f}uK".format(rms4),
    sub=[nrow, ncol, 2],
)
hp.graticule(**grargs)

hp.gnomview(
    (m4 - m1)[0] * 1e6,
    **args,
    title="difference, rms={:.1f}uK".format(rms14),
    sub=[nrow, ncol, 3],
)
hp.graticule(**grargs)

cl1 = hp.anafast(m1[0] * mask, lmax=1024)
cl4 = hp.anafast(m4[0] * mask, lmax=1024)
cld = hp.anafast((m4-m1)[0] * mask, lmax=1024)
ell = np.arange(1025)
ax = fig.add_subplot(nrow, ncol, 4)
ax.loglog(ell[1:], cl1[1:] * 1e12, label="thinfp=1")
ax.loglog(ell[1:], cl4[1:] * 1e12, label="thinfp=4")
ax.loglog(ell[1:], cld[1:] * 1e12, label="diff")
ax.set_xlabel(r"Multipole, $\ell$")
ax.set_ylabel(r"C$_\ell$ [$\mu$K$^2$]")
plt.legend(loc="best")

plt.savefig("atmosphere_vs_ndet.png")
