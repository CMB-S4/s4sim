import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from toast.pixels_io_healpix import read_healpix, write_healpix


telescope = "splat"
TELESCOPE = "LAT2_SPLAT"

band = "f150"
m1 = read_healpix(f"maps/single_obs_{telescope}_{band}.h5", None)
m2_1 = read_healpix(f"maps/single_pole1_{telescope}_{band}.h5", None)
m2_2 = read_healpix(f"maps/single_pole2_{telescope}_{band}.h5", None)
m2_3 = read_healpix(f"maps/single_pole3_{telescope}_{band}.h5", None)
m3 = read_healpix(f"/global/cfs/cdirs/cmbs4/dc/dc0/staging/noise_sim/outputs_rk/coadd/{telescope}/coadd_{telescope}_{band}_001of001_map.fits", None)

for m in [m1, m2_1, m2_2, m2_3, m3]:
    m *= 1e6
    m[m == 0] = hp.UNSEEN

args = {"rot" : [40, -55], "half_sky" : True, "unit" : "$\mu$K", "cmap" : "bwr", "xsize" : 800}
nrow, ncol = 1, 5
scale = 5 / 3 * .85 / 1.0
fig = plt.figure(figsize=[scale * 4 * ncol, scale * 5 * nrow])
amp = 200
hp.orthview(m1[0], min=-amp, max=amp, sub=[nrow, ncol, 1], title="Single Observation", **args)
hp.orthview(m2_1[0], min=-amp, max=amp, sub=[nrow, ncol, 2], title="Single Pass - Field 1", **args)
hp.orthview(m2_2[0], min=-amp, max=amp, sub=[nrow, ncol, 3], title="Single Pass - Field 2", **args)
hp.orthview(m2_3[0], min=-amp, max=amp, sub=[nrow, ncol, 4], title="Single Pass - Field 3", **args)
amp = 2
hp.orthview(m3[0], min=-amp, max=amp, sub=[nrow, ncol, 5], title="10 years", **args)

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

nrow, ncol = 1, 3
fig = plt.figure(figsize=[4 * ncol, 5 * nrow])
hp.orthview(mask, sub=[nrow, ncol, 1], title="Mask", **args)
amp = 300
hp.orthview(cmb, min=-amp, max=amp, sub=[nrow, ncol, 2], title="CMB", **args)
hp.orthview(fg, min=-amp, max=amp, sub=[nrow, ncol, 3], title="Foregrounds", **args)

fig.savefig(f"mask_and_signal_{telescope}_{band}.png")
fig.savefig(f"mask_and_signal_{telescope}_{band}.pdf")
