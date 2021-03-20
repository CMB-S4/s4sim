# Adapted from Joel Myers' notebook for CMB-S4 program level requirements

import sys
import pickle
import numpy as np
import pylab
from pylab import *
import os
import matplotlib.cm as cm

# pylab.rc('text', usetex=True)
# pylab.rc('font', family='serif')

deconvolve_beam = True

with open("fiducialSpectra.pkl", "rb") as f:
    cambPowersFid = pickle.load(f, encoding="latin1")

fiducial_ell = cambPowersFid["lensed"]["l"]
fiducial_TT = cambPowersFid["lensed"]["dl_TT"]
fiducial_EE = cambPowersFid["lensed"]["dl_EE"]
fiducial_BB = cambPowersFid["lensed"]["dl_BB"]


def get_bl(fwhm_arcmin, ells):
    fwhm_radians = np.radians(fwhm_arcmin / 60)
    sigma = fwhm_radians / np.sqrt(8 * np.log(2))
    bl = np.exp(ells * (ells + 1) * sigma ** 2)
    return bl


def get_nl(
    noiseval, el, beamval=None, use_beam_window=0, uk_to_K=0, elknee_t=-1, alpha_knee=0
):

    if uk_to_K:
        noiseval = noiseval / 1e6

    if use_beam_window:
        bl = get_bl(beamval, el)

    delta_T_radians = noiseval * np.radians(1 / 60)
    nl = np.tile(delta_T_radians ** 2, int(max(el)) + 1)

    nl = np.asarray([nl[int(l)] for l in el])

    if use_beam_window:
        nl *= bl

    if elknee_t != -1:
        nl = np.copy(nl) * (1 + (elknee_t * 1 / el) ** alpha_knee)

    return nl

band2freq = {
    "ULF" : 20,
    "LFL1" : 27,
    "LFL2" : 39,
    "MFL1" : 93,
    "MFL2" : 145,
    "HFL1" : 225,
    "HFL2" : 278,
    "LFS1" : 30,
    "LFS2" : 40,
    "MFLS1" : 85,
    "MFLS2" : 145,
    "MFHS1" : 95,
    "MFHS2" : 155,
    "HFS1" : 220,
    "HFS2" : 270,
}

Chile_LAT = {
    # freq: [beam_arcmins, white_noise_T, elknee_T, alphaknee_T, whitenoise_P, elknee_P, alphaknee_P]
    27: [7.4, 21.34, 415.0, 3.5, 30.23, 700, 1.4],
    39: [5.1, 11.67, 391.0, 3.5, 16.53, 700, 1.4],
    93: [2.2, 1.89, 1932.0, 3.5, 2.68, 700, 1.4],
    145: [1.4, 2.09, 3917.0, 3.5, 2.96, 700, 1.4],
    225: [1.0, 6.90, 6740.0, 3.5, 9.78, 700, 1.4],
    278: [0.9, 16.88, 6792.0, 3.5, 23.93, 700, 1.4],
}

Pole_LAT = {
    # freq: [beam_arcmins, white_noise_T, elknee_T, alphaknee_T, whitenoise_P, elknee_P, alphaknee_P]
    20: [11.0, 9.31, 1200, 4.2, 13.16, 150, 2.7],
    27: [8.4, 4.6, 1200.0, 4.2, 6.5, 150, 2.7],
    39: [5.8, 2.94, 1200.0, 4.2, 4.15, 150, 2.7],
    93: [2.5, 0.45, 1200.0, 4.2, 0.63, 150, 2.6],
    145: [1.6, 0.41, 1900.0, 4.1, 0.59, 200, 2.2],
    225: [1.1, 1.29, 2100.0, 4.1, 1.83, 200, 2.2],
    278: [1.0, 3.07, 2100.0, 3.9, 4.34, 200, 2.2],
}

Pole_SAT = {
    # freq: [beam_arcmins, white_noise_T, elknee_T, alphaknee_T, whitenoise_E, elknee_E, alphaknee_E, whitenoise_B, elknee_B, alphaknee_B]
    30: [72.8, 5.64, 150, 4.4, 3.74, 60, 2.2, 3.53, 60, 1.7],
    40: [72.8, 7.14, 150.0, 4.4, 4.73, 60, 2.2, 4.46, 60, 1.7],
    85: [25.5, 1.41, 150.0, 4.4, 0.93, 60, 2.2, 0.88, 60, 1.7],
    95: [22.7, 1.24, 150.0, 4.4, 0.82, 60, 2.2, 0.78, 60, 1.7],
    145: [25.5, 2.71, 230.0, 3.8, 1.25, 65, 3.1, 1.23, 60, 3.0],
    155: [22.7, 2.90, 230.0, 3.8, 1.34, 65, 3.1, 1.34, 60, 3.0],
    220: [13.0, 7.50, 230.0, 3.8, 3.48, 65, 3.1, 3.48, 60, 3.0],
    270: [13.0, 12.85, 230.0, 3.8, 8.08, 65, 3.1, 5.97, 60, 3.0],
}

lmax = 6000
ells = np.arange(2, lmax)
ellnorm = ells * (ells + 1) / (2 * np.pi)
NlTT_Chile_LAT = dict()
NlEE_Chile_LAT = dict()
NlTT_Pole_LAT = dict()
NlEE_Pole_LAT = dict()
NlTT_Pole_SAT = dict()
NlEE_Pole_SAT = dict()
NlBB_Pole_SAT = dict()
for freq in Chile_LAT.keys():
    NlTT_Chile_LAT[freq] = get_nl(
        noiseval=Chile_LAT[freq][1],
        el=ells,
        beamval=Chile_LAT[freq][0],
        use_beam_window=deconvolve_beam,
        elknee_t=Chile_LAT[freq][2],
        alpha_knee=Chile_LAT[freq][3],
    ) * ellnorm
    NlEE_Chile_LAT[freq] = get_nl(
        noiseval=Chile_LAT[freq][4],
        el=ells,
        beamval=Chile_LAT[freq][0],
        use_beam_window=deconvolve_beam,
        elknee_t=Chile_LAT[freq][5],
        alpha_knee=Chile_LAT[freq][6],
    ) * ellnorm


for freq in Pole_LAT.keys():
    NlTT_Pole_LAT[freq] = get_nl(
        noiseval=Pole_LAT[freq][1],
        el=ells,
        beamval=Pole_LAT[freq][0],
        use_beam_window=deconvolve_beam,
        elknee_t=Pole_LAT[freq][2],
        alpha_knee=Pole_LAT[freq][3],
    ) * ellnorm
    NlEE_Pole_LAT[freq] = get_nl(
        noiseval=Pole_LAT[freq][4],
        el=ells,
        beamval=Pole_LAT[freq][0],
        use_beam_window=deconvolve_beam,
        elknee_t=Pole_LAT[freq][5],
        alpha_knee=Pole_LAT[freq][6],
    ) * ellnorm

for freq in Pole_SAT.keys():
    NlTT_Pole_SAT[freq] = get_nl(
        noiseval=Pole_SAT[freq][1],
        el=ells,
        beamval=Pole_SAT[freq][0],
        use_beam_window=deconvolve_beam,
        elknee_t=Pole_SAT[freq][2],
        alpha_knee=Pole_SAT[freq][3],
    ) * ellnorm
    NlEE_Pole_SAT[freq] = get_nl(
        noiseval=Pole_SAT[freq][4],
        el=ells,
        beamval=Pole_SAT[freq][0],
        use_beam_window=deconvolve_beam,
        elknee_t=Pole_SAT[freq][5],
        alpha_knee=Pole_SAT[freq][6],
    ) * ellnorm

    NlBB_Pole_SAT[freq] = get_nl(
        noiseval=Pole_SAT[freq][7],
        el=ells,
        beamval=Pole_SAT[freq][0],
        use_beam_window=deconvolve_beam,
        elknee_t=Pole_SAT[freq][8],
        alpha_knee=Pole_SAT[freq][9],
    ) * ellnorm
