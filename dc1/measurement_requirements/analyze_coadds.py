import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
# from toast.pixels_io import read_healpix

import requirements as req


#multipanel = False
multipanel = True
#fileformat = "pdf"
fileformat = "png"

nside = 4096
lmax = 2 * nside
ell = np.arange(lmax + 1)

# These factors are now already included in the coadded maps
ntele = 1
nseason = 1

rootdir = "/global/cfs/cdirs/cmbs4/dc/dc1/staging/noise_sim/outputs_rk"
#rootdir = "outputs"

#for band in 90, 150:
for band in 30, 40, 90, 150, 220, 280,:
    tele = "chlat"
    TELE = "LAT0_CHLAT"

    fname_map = f"{rootdir}/coadd/{TELE}/coadd_{TELE}_f{band:03}_001of001_map.fits"
    fname_cov = f"{rootdir}/coadd/{TELE}/coadd_{TELE}_f{band:03}_001of001_cov.fits"
    fname_cl = f"outputs/cl/coadd_{TELE}_f{band:03}_001of001_cl.fits"

    fname_cl_cmb_in = f"../cmb_sim/outputs/cl/{TELE}/input_{TELE}_f{band:03}_cl.fits"
    fname_cl_cmb_out = f"../cmb_sim/outputs/cl/{TELE}/coadd_{TELE}_f{band:03}_cl.fits"
    print(f"Loading {fname_cl_cmb_in}", flush=True)
    cl_in = hp.read_cl(fname_cl_cmb_in)
    print(f"Loading {fname_cl_cmb_out}", flush=True)
    cl_out = hp.read_cl(fname_cl_cmb_out)
    tf = cl_out[:3] / cl_in[:3]
    tf[:, :2] = 1
    tf[:, 2000:] = 1

    if os.path.isfile(fname_cl):
        print(f"Loading {fname_cl}", flush=True)
        cl = hp.read_cl(fname_cl)
    else:
        outroot = os.path.dirname(fname_cl)
        os.makedirs(outroot, exist_ok=True)
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
        print(f"Writing {fname_cl}", flush=True)
        hp.write_cl(fname_cl, cl, overwrite=True)

    freq = {30 : 27, 40 : 39, 90 : 93, 150 : 145, 220 : 225, 280 : 278}[band]
    fwhm = req.Chile_LAT[freq][0]
    bl = req.get_bl(fwhm, ell)
    nltt = req.NlTT_Chile_LAT[freq]
    nlee = req.NlEE_Chile_LAT[freq]

    if multipanel:
        nrow, ncol = 1, 3
    else:
        nrow, ncol = 1, 1
    fig = plt.figure(figsize=[4 * ncol, 3 * nrow])
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
    ax.loglog(ell[2:], (ellnorm * cl[0] * scale)[2:], label=f"DC0")
    ax.loglog(ell[2:], (ellnorm * cl[0] * bl * scale / tf[0])[2:], label=f"TF-corrected DC0")
    ax.legend(loc="best")
    if not multipanel:
        fname_out = f"cl_comparison_noise_{band:03}_TT.{fileformat}"
        fig.savefig(fname_out)
        print(f"Wrote {fname_out}", flush=True)
        plt.clf()
        iplot = 0

    iplot += 1
    ax = fig.add_subplot(nrow, ncol, iplot)
    ax.set_title(f"EE {band}GHz")
    ax.set_xlabel("Multipole, $\ell$")
    ax.set_ylabel("D$\ell$ [$\mu$K$^2$]")
    ax.loglog(req.fiducial_ell, req.fiducial_EE, "k", label="CMB")
    ax.loglog(req.ells, nlee, label="requirement")
    ax.loglog(ell[2:], (ellnorm * cl[1] * scale)[2:], label=f"DC0")
    ax.loglog(ell[2:], (ellnorm * cl[1] * bl * scale / tf[1])[2:], label=f"DC0")
    if not multipanel:
        fname_out = f"cl_comparison_noise_{band:03}_EE.{fileformat}"
        fig.savefig(fname_out)
        print(f"Wrote {fname_out}", flush=True)
        plt.clf()
        iplot = 0

    iplot += 1
    ax = fig.add_subplot(nrow, ncol, iplot)
    ax.set_title(f"BB {band}GHz")
    ax.set_xlabel("Multipole, $\ell$")
    ax.set_ylabel("D$\ell$ [$\mu$K$^2$]")
    ax.loglog(req.fiducial_ell, req.fiducial_BB, "k", label="CMB")
    ax.loglog(req.ells, nlee, label="requirement")
    ax.loglog(ell[2:], (ellnorm * cl[2] * scale)[2:], label=f"DC0")
    ax.loglog(ell[2:], (ellnorm * cl[2] * bl * scale / tf[2])[2:], label=f"DC0")
    fig.tight_layout()
    if multipanel: 
        fname_out = f"cl_comparison_noise_{band:03}.{fileformat}"
        fig.savefig(fname_out)
        print(f"Wrote {fname_out}", flush=True)
    else:
        fname_out = f"cl_comparison_noise_{band:03}_BB.{fileformat}"
        fig.savefig(fname_out)
        print(f"Wrote {fname_out}", flush=True)
