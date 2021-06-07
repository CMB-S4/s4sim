# This script takes the simulated hitmaps and produces boolean coverage maps in galactic coordinates

import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np


hmap_in = hp.read_map("out-old/00000000/chile_noise_LAT_MFL1_filtered_telescope_all_time_all_hmap.fits")

nside_in = hp.get_nside(hmap_in)
npix_in = 12 * nside_in ** 2
nside_out = 256
npix_out = 12 * nside_out ** 2

hmap_out = np.zeros(npix_out, dtype=bool)

rot = hp.rotator.Rotator(coord="CG")

pix_in = np.arange(npix_in, dtype=np.int32)[hmap_in != 0]
vec_in = hp.pix2vec(nside_in, pix_in)
vec_out = rot(vec_in)
pix_out = hp.vec2pix(nside_out, *vec_out)
hmap_out[pix_out] = True

hp.write_map("CMBS4_Chile_coverage.fits", hmap_out, coord="G", overwrite=True)

#######


hmap_in = hp.read_map("out-old/00000000/pole_noise_SAT_MFHS1_telescope_all_time_all_hits.fits.gz")

nside_in = hp.get_nside(hmap_in)
npix_in = 12 * nside_in ** 2
nside_out = 256
npix_out = 12 * nside_out ** 2

hmap_out = np.zeros(npix_out, dtype=bool)

rot = hp.rotator.Rotator(coord="CG")

pix_in = np.arange(npix_in, dtype=np.int32)[hmap_in != 0]
vec_in = hp.pix2vec(nside_in, pix_in)
vec_out = rot(vec_in)
pix_out = hp.vec2pix(nside_out, *vec_out)
hmap_out[pix_out] = True

hp.write_map("CMBS4_Pole_coverage.fits", hmap_out, coord="G", overwrite=True)
