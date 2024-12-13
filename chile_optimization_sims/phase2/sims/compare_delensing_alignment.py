import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np


sat = hp.read_map("scaled_outputs/sun90_f155_90years_depth.fits")
delensing = hp.read_map("scaled_outputs/lat_delensing_f150_36years_depth.fits")
core = hp.read_map("scaled_outputs/lat_delensing_core_f150_36years_depth.fits")

# Normalize

for m in [sat, delensing, core]:
    good = m != 0
    m[good] /= np.amin(m[good])
    m[m == 0] = hp.UNSEEN

# Plot
nrow, ncol = 1, 3
fig = plt.figure(figsize=[ncol * 6, nrow * 4])

args = {
    "xsize" : 1600,
    "min" : 1,
    "max" : 2,
    "cmap" : "magma",
    "unit" : "Relative depth",
}
hp.mollview(sat, **args, title="SAT", sub=[nrow, ncol, 1])
hp.mollview(delensing, **args, title="delensing", sub=[nrow, ncol, 2])
hp.mollview(core, **args, title="delensing_core", sub=[nrow, ncol, 3])

nrow, ncol = 1, 3
fig = plt.figure(figsize=[ncol * 4, nrow * 4])

args.update({
    "xsize" : 800,
    "rot" : [45, -45],
    "reso" : 5.0,
})
gr = 10
grcolor = "white"
hp.gnomview(sat, **args, title="SAT", sub=[nrow, ncol, 1])
hp.graticule(gr, color=grcolor)
hp.gnomview(delensing, **args, title="delensing", sub=[nrow, ncol, 2])
hp.graticule(gr, color=grcolor)
hp.gnomview(core, **args, title="delensing_core", sub=[nrow, ncol, 3])
hp.graticule(gr, color=grcolor)
