import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from toast.pixels_io import read_healpix

import requirements as req


nside = 4096
lmax = 2 * nside
ell = np.arange(lmax + 1)

ntele = 2
nseason = 7

tele = "chlat"

for band in 90, 150:
    fname_input_map = f"input_maps/cmb.{tele}.f{band:03}.h5"
    fname_map = f"outputs/coadd/{tele}/coadd_{tele}_f{band:03}_map.fits"
    fname_cov = f"outputs/coadd/{tele}/coadd_{tele}_f{band:03}_cov.fits"
    fname_cl = f"outputs/coadd/{tele}/coadd_{tele}_f{band:03}_cl.fits"
    fname_input_cl = f"outputs/coadd/{tele}/input_{tele}_f{band:03}_cl.fits"

    if os.path.isfile(fname_cl):
        cl = hp.read_cl(fname_cl)
        cl_in = hp.read_cl(fname_input_cl)
    else:
        print(f"Loading {fname_map}")
        m = hp.read_map(fname_map, None)
        print(f"Loading {fname_cov}")
        w = hp.read_map(fname_cov)
        # Discard 1% of the noisiest pixels
        good = w > 0
        sorted_w = np.sort(w[good])
        ngood = sorted_w.size
        lim = sorted_w[int(.99 * ngood)]
        mask = np.logical_and(w > 0, w < lim)
        print("Measuring C_ell")
        m[0] = hp.remove_dipole(m[0] * mask, bad=0)
        cl = hp.anafast(m * mask, lmax=lmax, iter=0)
        fsky = np.sum(mask) / mask.size
        cl /= fsky
        print(f"Writing {fname_cl}")
        hp.write_cl(fname_cl, cl, overwrite=True)

        print(f"Loading {fname_input_map}")
        inmap = read_healpix(fname_input_map, None)
        inmap[0] = hp.remove_dipole(inmap[0] * mask, bad=0)
        print("Measuring input C_ell")
        cl_in = hp.anafast(inmap * mask, lmax=lmax, iter=0)
        cl_in /= fsky
        print(f"Writing {fname_cl}")
        hp.write_cl(fname_input_cl, cl_in, overwrite=True)

    freq = {30 : 27, 40 : 39, 90 : 93, 150 : 145, 220 : 225, 280 : 278}[band]
    fwhm = req.Chile_LAT[freq][0]
    bl = req.get_bl(fwhm, ell)
    nltt = req.NlTT_Chile_LAT[freq]
    nlee = req.NlEE_Chile_LAT[freq]

    nrow, ncol = 2, 2
    fig = plt.figure(figsize=[6 * ncol, 6 * nrow])
    ell = np.arange(lmax + 1)
    ellnorm = ell * (ell + 1) / (2 * np.pi) * 1e12
    scale = 1 # / ntele / nseason

    iplot = 1
    ax = fig.add_subplot(nrow, ncol, iplot)
    ax.set_title(f"TT {band} / {freq}GHz")
    ax.set_xlabel("Multipole, $\ell$")
    ax.set_ylabel("D$\ell$ [$\mu$K$^2$]")
    ax.loglog(req.fiducial_ell, req.fiducial_TT, "k", label="CMB")
    ax.loglog(req.ells, nltt, label="requirement")
    ax.loglog(ell[2:], (ellnorm * cl_in[0] * bl * scale)[2:], label=f"DC1 input")
    ax.loglog(ell[2:], (ellnorm * cl[0] * bl * scale)[2:], label=f"DC1 output")

    iplot += 1
    ax = fig.add_subplot(nrow, ncol, iplot)
    ax.set_title(f"EE {band} / {freq}GHz")
    ax.set_xlabel("Multipole, $\ell$")
    ax.set_ylabel("D$\ell$ [$\mu$K$^2$]")
    ax.loglog(req.fiducial_ell, req.fiducial_EE, "k", label="CMB")
    ax.loglog(req.ells, nlee, label="requirement")
    ax.loglog(ell[2:], (ellnorm * cl_in[1] * bl * scale)[2:], label=f"DC1 input")
    ax.loglog(ell[2:], (ellnorm * cl[1] * bl * scale)[2:], label=f"DC1 output")

    iplot += 1
    ax = fig.add_subplot(nrow, ncol, iplot)
    ax.set_title(f"BB {band} / {freq}GHz")
    ax.set_xlabel("Multipole, $\ell$")
    ax.set_ylabel("D$\ell$ [$\mu$K$^2$]")
    ax.loglog(req.fiducial_ell, req.fiducial_BB, "k", label="CMB")
    ax.loglog(req.ells, nlee, label="requirement")
    ax.loglog(ell[2:], (ellnorm * cl_in[2] * bl * scale)[2:], label=f"DC1 input")
    ax.loglog(ell[2:], (ellnorm * cl[2] * bl * scale)[2:], label=f"DC1 output")

    fname_out = f"cl_comparison_cmb_{band:03}.png"
    fig.savefig(fname_out)
    print(f"Wrote {fname_out}")
