import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from toast.pixels_io_healpix import read_healpix, write_healpix


telescope = "chlat"
TELESCOPE = "LAT0_CHLAT"

band = "f150"
m1 = read_healpix(f"maps/single_obs_{telescope}_{band}.h5", None)
m2 = read_healpix(f"maps/single_day_{telescope}_{band}.h5", None)
m3 = read_healpix(f"/global/cfs/cdirs/cmbs4/dc/dc0/staging/noise_sim/outputs_rk/coadd/{TELESCOPE}/coadd_{TELESCOPE}_{band}_001of001_map.fits", None)

for m in [m1, m2, m3]:
    m *= 1e6
    m[m == 0] = hp.UNSEEN

args = {"unit" : "$\mu$K", "cmap" : "bwr", "xsize" : 1600}
nrow, ncol = 1, 3
fig = plt.figure(figsize=[4 * ncol, 3 * nrow])
amp = 1000
hp.mollview(m1[0], min=-amp, max=amp, sub=[nrow, ncol, 1], title="Single Observation", **args)
hp.mollview(m2[0], min=-amp, max=amp, sub=[nrow, ncol, 2], title="Single Day", **args)
amp = 10
hp.mollview(m3[0], min=-amp, max=amp, sub=[nrow, ncol, 3], title="7 years", **args)

fig.savefig(f"obs_maps_{telescope}_{band}.png")
fig.savefig(f"obs_maps_{telescope}_{band}.pdf")

mask = read_healpix(f"../foreground_sim/input_maps/mask_mediumcomplexity.{telescope}.{band}.h5", None, dtype=float)[0]
cmb = read_healpix(f"../cmb_sim/input_maps/cmb.{telescope}.{band}.h5", None)[0]
cmb = hp.remove_dipole(cmb)
fg = read_healpix(f"../foreground_sim/input_maps/foreground_mediumcomplexity.{telescope}.{band}.h5", None)[0]

mask[m3[0] == hp.UNSEEN] = hp.UNSEEN
for m in [cmb, fg]:
    m *= 1e6
    m[m3[0] == hp.UNSEEN] = hp.UNSEEN

fig = plt.figure(figsize=[4 * ncol, 3 * nrow])
hp.mollview(mask, sub=[nrow, ncol, 1], title="Mask", **args)
amp = 300
hp.mollview(cmb, min=-amp, max=amp, sub=[nrow, ncol, 2], title="CMB", **args)
hp.mollview(fg, min=-amp, max=amp, sub=[nrow, ncol, 3], title="Foregrounds", **args)

fig.savefig(f"mask_and_signal_{telescope}_{band}.png")
fig.savefig(f"mask_and_signal_{telescope}_{band}.pdf")
