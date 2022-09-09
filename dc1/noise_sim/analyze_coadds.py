import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
# from toast.pixels_io import read_healpix

import requirements as req

multipanel = False

nside = 4096
lmax = 2 * nside
ell = np.arange(lmax + 1)

ntele = 2
nseason = 7

#for band in 90, 150:
for band in 30, 40:
    tele = "chlat"
    TELE = "LAT0_CHLAT"

    fname_map = f"outputs/coadd/{TELE}/coadd_{TELE}_f{band:03}_map.fits"
    fname_cov = f"outputs/coadd/{TELE}/coadd_{TELE}_f{band:03}_cov.fits"
    fname_cl = f"outputs/coadd/{TELE}/coadd_{TELE}_f{band:03}_cl.fits"

    fname_cl_cmb_in = f"../cmb_sim/outputs/coadd/{TELE}/input_{TELE}_f{band:03}_cl.fits"
    fname_cl_cmb_out = f"../cmb_sim/outputs/coadd/{TELE}/coadd_{TELE}_f{band:03}_cl.fits"
    cl_in = hp.read_cl(fname_cl_cmb_in)
    cl_out = hp.read_cl(fname_cl_cmb_out)
    tf = cl_out[:3] / cl_in[:3]
    tf[:, :2] = 1
    tf[:, 2000:] = 1

    if os.path.isfile(fname_cl):
        cl = hp.read_cl(fname_cl)
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

    freq = {30 : 27, 40 : 39, 90 : 93, 150 : 145, 220 : 225, 280 : 278}[band]
    fwhm = req.Chile_LAT[freq][0]
    bl = req.get_bl(fwhm, ell)
    nltt = req.NlTT_Chile_LAT[freq]
    nlee = req.NlEE_Chile_LAT[freq]

    if multipanel:
        nrow, ncol = 2, 2
    else:
        nrow, ncol = 1, 1
    fig = plt.figure(figsize=[6 * ncol, 6 * nrow])
    ell = np.arange(lmax + 1)
    ellnorm = ell * (ell + 1) / (2 * np.pi) * 1e12
    scale = 1 / ntele / nseason

    iplot = 1
    ax = fig.add_subplot(nrow, ncol, iplot)
    ax.set_title(f"TT {band}GHz")
    ax.set_xlabel("Multipole, $\ell$")
    ax.set_ylabel("D$\ell$ [$\mu$K$^2$]")
    ax.loglog(req.fiducial_ell, req.fiducial_TT, "k", label="fiducial")
    ax.loglog(req.ells, nltt, label="requirement")
    ax.loglog(ell[2:], (ellnorm * cl[0] * scale)[2:], label=f"DC1")
    ax.loglog(ell[2:], (ellnorm * cl[0] * bl * scale / tf[0])[2:], label=f"TF-corrected DC1")
    ax.legend(loc="best")
    if not multipanel:
        fname_out = f"cl_comparison_noise_{band:03}_TT.png"
        fig.savefig(fname_out)
        print(f"Wrote {fname_out}")
        plt.clf()
        iplot = 0

    iplot += 1
    ax = fig.add_subplot(nrow, ncol, iplot)
    ax.set_title(f"EE {band}GHz")
    ax.set_xlabel("Multipole, $\ell$")
    ax.set_ylabel("D$\ell$ [$\mu$K$^2$]")
    ax.loglog(req.fiducial_ell, req.fiducial_EE, "k", label="CMB")
    ax.loglog(req.ells, nlee, label="requirement")
    ax.loglog(ell[2:], (ellnorm * cl[1] * scale)[2:], label=f"DC1")
    ax.loglog(ell[2:], (ellnorm * cl[1] * bl * scale / tf[1])[2:], label=f"DC1")
    if not multipanel:
        fname_out = f"cl_comparison_noise_{band:03}_EE.png"
        fig.savefig(fname_out)
        print(f"Wrote {fname_out}")
        plt.clf()
        iplot = 0

    iplot += 1
    ax = fig.add_subplot(nrow, ncol, iplot)
    ax.set_title(f"BB {band}GHz")
    ax.set_xlabel("Multipole, $\ell$")
    ax.set_ylabel("D$\ell$ [$\mu$K$^2$]")
    ax.loglog(req.fiducial_ell, req.fiducial_BB, "k", label="CMB")
    ax.loglog(req.ells, nlee, label="requirement")
    ax.loglog(ell[2:], (ellnorm * cl[2] * scale)[2:], label=f"DC1")
    ax.loglog(ell[2:], (ellnorm * cl[2] * bl * scale / tf[2])[2:], label=f"DC1")
    if multipanel: 
        fname_out = f"cl_comparison_noise_{band:03}.png"
        fig.savefig(fname_out)
        print(f"Wrote {fname_out}")
    else:
        fname_out = f"cl_comparison_noise_{band:03}_BB.png"
        fig.savefig(fname_out)
        print(f"Wrote {fname_out}")
