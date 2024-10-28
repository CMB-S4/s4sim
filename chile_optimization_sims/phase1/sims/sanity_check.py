import glob
import os
import pickle
import sys

import healpy as hp
import numpy as np


ntube = {
    "f030" : 3,
    "f040" : 3,
    "f085" : 9,
    "f095" : 9,
    "f145" : 9,
    "f155" : 9,
    "f220" : 6,
    "f280" : 6,
}

f_total = {
    "f030" : 0.17,
    "f040" : 0.17,
    "f085" : 0.17,
    "f095" : 0.17,
    "f145" : 0.16,
    "f155" : 0.16,
    "f220" : 0.10,
    "f280" : 0.07,
}

f_sensitivity = {
    "f030" : 0.8 * 0.8,
    "f040" : 0.8 * 0.8,
    "f085" : 0.8 * 0.8,
    "f095" : 0.8 * 0.8,
    "f145" : 0.8 * 0.65,
    "f155" : 0.8 * 0.65,
    "f220" : 0.8 * 0.65,
    "f280" : 0.8 * 0.65,
}

depths = {}

for freq in f_total:
    ftot = f_total[freq]
    fsen = f_sensitivity[freq]

    fname_fp = glob.glob(f"focalplane_*_SAT_{freq}_*.h5_thinfp=None_fsample=1.0.pck")
    if len(fname_fp) == 0:
        fname_fp = glob.glob(f"old_focalplanes/focalplane_*_SAT_{freq}_*.h5_thinfp=None_fsample=1.0.pck")
    fname_fp = fname_fp[0]
    # print(f"{freq}: {fname_fp}")
    fp = pickle.load(open(fname_fp, "rb"))
    ndet = len(fp.detectors)
    net = np.mean(fp.detector_data["psd_net"].to_value()) * 1e6
    anet = net / np.sqrt(ndet)

    # Relative depth number to compare
    depth = anet / np.sqrt(ftot * fsen * ntube[freq])
    depths[freq] = depth

    print(f"{freq:3} : {ftot:.3f} * {fsen:.3f} = {ftot * fsen:.3f}, array_net = {anet:6.3f}, depth = {depth:7.3f}")

freqs = sorted(f_total.keys())
for ifreq1, freq1 in enumerate(freqs):
    depth1 = depths[freq1]
    if True:
        fname1 = glob.glob(f"scaled_outputs/sat_{freq1}_*_depth.fits")[0]
        m1 = hp.read_map(fname1) * np.sqrt(2)  # T->P depth
    else:
        fname1 = f"outputs/sat/{freq1}/mapmaker_cov.fits"
        # fname1 = glob.glob(f"scaled_outputs/sat_{freq1}_*_cov.fits")[0]
        m1 = np.sqrt(hp.read_map(fname1) * 2) * 1e6  # T->P depth
    lim1 = np.sort(m1[m1 != 0])[int(m1.size * 0.03)]
    best1 = m1[np.logical_and(m1 != 0, m1 < lim1)]
    mean1 = np.mean(best1)
    for freq2 in freqs[ifreq1 + 1:]:
        depth2 = depths[freq2]
        if True:
            fname2 = glob.glob(f"scaled_outputs/sat_{freq2}_*_depth.fits")[0]
            m2 = hp.read_map(fname2) * np.sqrt(2)  # T->P depth
        else:
            fname2 = f"outputs/sat/{freq2}/mapmaker_cov.fits"
            # fname2 = glob.glob(f"scaled_outputs/sat_{freq2}_*_cov.fits")[0]
            m2 = np.sqrt(hp.read_map(fname2) * 2) * 1e6  # T->P depth
        lim2 = np.sort(m2[m2 != 0])[int(m2.size * 0.03)]
        best2 = m2[np.logical_and(m2 != 0, m2 < lim2)]
        mean2 = np.mean(best2)
        ratio1 = depth1 / depth2
        ratio2 = mean1 / mean2
        print(f"Depth ratio: {freq1} / {freq2} = {depth1:7.3f} / {depth2:7.3f} = {ratio1:.3f} vs. {mean1:6.3f} / {mean2:6.3f} = {ratio2:.3f} ({ratio1 / ratio2:.3f})")
