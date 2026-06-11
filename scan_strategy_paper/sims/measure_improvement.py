import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np


fname_classic = "outputs/sat_classic/SAT_f090/mapmaker_invcov.fits"
fname_deep = "outputs/sat_deep/SAT_f090/mapmaker_invcov.fits"

classic = hp.read_map(fname_classic)
deep = hp.read_map(fname_deep)

for pc in 0, 50, 75, 90, 95, 97, 98, 99:
    limit_classic = np.percentile(classic, pc)
    limit_deep = np.percentile(deep, pc)

    mask_classic = classic > limit_classic
    mask_deep = deep > limit_deep

    weight_classic = np.sum(classic[mask_classic])
    weight_deep = np.sum(deep[mask_deep])

    print(f"Survey weight ratio (best {100 - pc:3}%) deep / classic = {weight_deep / weight_classic:.3f}")
