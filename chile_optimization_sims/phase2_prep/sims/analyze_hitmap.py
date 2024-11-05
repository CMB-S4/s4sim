import os
import sys

import healpy as hp
import numpy as np


if len(sys.argv) == 1:
    print(f"Usage: {sys.argv[0]} <hitmap>")
    sys.exit()

fname = sys.argv[1]
if not os.path.isfile(fname):
    print(f"No such hitmap {fname}")
    sys.exit()

hits = hp.read_map(fname)
npix = hits.size
good = hits > 0
dhit = hits[good].astype(float)

if "invcov" in fname:
    print(f"Interpreting map as inverse covariance")
elif "cov" in fname:
    print(f"Interpreting map as covariance")
    dhit = 1 / dhit
elif "depth" in fname:
    print(f"Interpreting map as noise depth")
    dhit = 1 / dhit**2

moment0 = np.sum(good) / npix
moment1 = np.sum(dhit) / npix
moment2 = np.sum(dhit ** 2) / npix
moment4 = np.sum(dhit ** 4) / npix
fraw = moment0
fnoise = moment1 ** 2 / moment2
fsignal = moment2 ** 2 / moment4

print(f"{fname}")
print(f"             Raw fsky: {fraw:.3f}")
print(f" Noise-dominated fsky: {fnoise:.3f}")
print(f"Signal-dominated fsky: {fsignal:.3f}")
