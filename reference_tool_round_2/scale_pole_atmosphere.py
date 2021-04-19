# Change the atmospheric gain in the simulated Pole maps to half of the original gain

from glob import glob
import os
import sys

import astropy.io.fits as pf
import healpy as hp
import numpy as np


scale = 0.5

fnames = sorted(glob("out/*/pole_atmosphere_*fits"))
for fname in fnames:
    print(f"Loading {fname}")
    hdulist = pf.open(fname, "update", memmap=False)
    hdu = hdulist[1]
    if "atmgain" in hdu.header:
        print(f"{fname} is already scaled")
        hdulist.close()
        continue
    good = hdu.data.field(0) != hp.UNSEEN
    for i in range(3):
        hdu.data.field(i)[good] *= scale
    hdu.header["atmgain"] = scale
    hdulist.close()
    print(f"Wrote {fname}")
