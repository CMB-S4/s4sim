import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np


d0 = hp.read_map("outputs/lat_wide/f090/3/mapmaker_cov.fits")
d1 = hp.read_map("outputs/lat_wide/f090/8/mapmaker_cov.fits")

good = d0 != 0
ratio = np.zeros_like(d0)
ratio[good] = np.sqrt(d1[good] / d0[good])
ratio[ratio==0] = hp.UNSEEN
med = np.median(ratio[good])
hp.mollview(ratio, title=f"90GHz, depth(8mm) / depth(3mm), median = {med:.3f}", unit="depth ratio", xsize=1600)

plt.savefig("pwv_depth_ratio_90GHz.png")


d0 = hp.read_map("outputs/lat_wide/f220/2/mapmaker_cov.fits")
d1 = hp.read_map("outputs/lat_wide/f220/8/mapmaker_cov.fits")

good = d0 != 0
ratio = np.zeros_like(d0)
ratio[good] = np.sqrt(d1[good] / d0[good])
ratio[ratio==0] = hp.UNSEEN
med = np.median(ratio[good])
hp.mollview(ratio, title=f"220GHz, depth(8mm) / depth(2mm), median = {med:.3f}", unit="depth ratio", xsize=1600)

plt.savefig("pwv_depth_ratio_220GHz.png")
