import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from toast.pixels_io_healpix import read_healpix

import requirements as req


multipanel = True

nside = 4096
lmax = 2 * nside
ell = np.arange(lmax + 1)

rootdir = "/global/cfs/cdirs/cmbs4/dc/dc0/staging"

for band in 30, 40, 90, 150, 220, 280:
#for band in 30, 40:
    tele = "chlat"
    TELE = "LAT0_CHLAT"

    fname_input = f"../cmb_sim/input_maps/cmb.chlat.f{band:03}.h5"
    fname_map = f"{rootdir}/multimap_sim/outputs_rk/coadd/{TELE}/coadd_{TELE}_f{band:03}_unlensed_cmb_001of001_map.fits"
    fname_mask = f"/global/cfs/cdirs/cmbs4/dc/dc0/masks/mask_mediumcomplexity.{tele}.f{band:03}.fits"
    fname_cov = f"{rootdir}/noise_sim/outputs_rk/coadd/{TELE}/coadd_{TELE}_f{band:03}_001of001_cov.fits"

    fname_cl_in = f"outputs/cl/{TELE}/input_cmb_{TELE}_f{band:03}_cl.fits"
    fname_cl_out = f"outputs/cl/{TELE}/coadd_cmb_{TELE}_f{band:03}_cl.fits"

    if not os.path.isfile(fname_cl_in) or not os.path.isfile(fname_cl_out):
        for path in fname_cl_in, fname_cl_out:
            outroot = os.path.dirname(path)
            os.makedirs(outroot, exist_ok=True)
        print(f"Loading {fname_map}", flush=True)
        input_map = read_healpix(fname_input)
        output_map = hp.read_map(fname_map, None)

        # Read and invert the processing mask
        print(f"Loading {fname_mask}")
        mask = np.logical_not(hp.read_map(fname_mask, dtype=bool))

        # Discard 1% of the noisiest pixels
        print(f"Loading {fname_cov}", flush=True)
        w = hp.read_map(fname_cov)
        good = w > 0
        sorted_w = np.sort(w[good])
        ngood = sorted_w.size
        lim = sorted_w[int(.99 * ngood)]
        mask[w == 0] = False
        mask[w > lim] = False

        print("Measuring C_ell", flush=True)
        input_map[0] = hp.remove_dipole(input_map[0] * mask, bad=0)
        output_map[0] = hp.remove_dipole(output_map[0] * mask, bad=0)
        cl_in = hp.anafast(input_map * mask, lmax=lmax, iter=0)
        cl_out = hp.anafast(output_map * mask, lmax=lmax, iter=0)
        fsky = np.sum(mask) / mask.size
        cl_in /= fsky
        cl_out /= fsky
        print(f"Writing {fname_cl_in}", flush=True)
        hp.write_cl(fname_cl_in, cl_in, overwrite=True)
        print(f"Writing {fname_cl_out}", flush=True)
        hp.write_cl(fname_cl_out, cl_out, overwrite=True)
    else:
        print(f"Loading {fname_cl_in}", flush=True)
        cl_in = hp.read_cl(fname_cl_in)
        print(f"Loading {fname_cl_out}", flush=True)
        cl_out = hp.read_cl(fname_cl_out)

    freq = {30 : 27, 40 : 39, 90 : 93, 150 : 145, 220 : 225, 280 : 278}[band]
    fwhm = req.Chile_LAT[freq][0]
    bl = req.get_bl(fwhm, ell)
    bl = 1 / bl[:req.fiducial_TT.size]
    nltt = req.NlTT_Chile_LAT[freq]
    nlee = req.NlEE_Chile_LAT[freq]

    if multipanel:
        nrow, ncol = 2, 2
    else:
        nrow, ncol = 1, 1
    fig = plt.figure(figsize=[6 * ncol, 6 * nrow])
    ell = np.arange(lmax + 1)
    ellnorm = ell * (ell + 1) / (2 * np.pi) * 1e12
    scale = 1 # / ntele / nseason

    iplot = 1
    ax = fig.add_subplot(nrow, ncol, iplot)
    ax.set_title(f"TT {band}GHz")
    ax.set_xlabel("Multipole, $\ell$")
    ax.set_ylabel("D$\ell$ [$\mu$K$^2$]")
    ax.loglog(req.fiducial_ell, req.fiducial_TT, "k", label="fiducial")
    ax.loglog(req.fiducial_ell, req.fiducial_TT * bl, "k--", label="fiducial x B$_\ell$")
    ax.loglog(req.ells, nltt, label="requirement")
    ax.loglog(ell[2:], (ellnorm * cl_in[0])[2:], label=f"C$_\ell$ in")
    ax.loglog(ell[2:], (ellnorm * cl_out[0])[2:], label=f"C$_\ell$ out")
    ax.legend(loc="best")
    if not multipanel:
        fname_out = f"cl_comparison_CMB_{band:03}_TT.png"
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
    ax.loglog(req.fiducial_ell, req.fiducial_EE * bl, "k--", label="fiducial x B$_\ell$")
    ax.loglog(req.ells, nlee, label="requirement")
    ax.loglog(ell[2:], (ellnorm * cl_in[1])[2:], label=f"C$_\ell$ in")
    ax.loglog(ell[2:], (ellnorm * cl_out[1])[2:], label=f"C$_\ell$ out")
    if not multipanel:
        fname_out = f"cl_comparison_CMB_{band:03}_EE.png"
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
    ax.loglog(req.fiducial_ell, req.fiducial_BB * bl, "k--", label="fiducial x B$_\ell$")
    ax.loglog(req.ells, nlee, label="requirement")
    ax.loglog(ell[2:], (ellnorm * cl_in[2])[2:], label=f"C$_\ell$ in")
    ax.loglog(ell[2:], (ellnorm * cl_out[2])[2:], label=f"C$_\ell$ out")
    if multipanel: 
        fname_out = f"cl_comparison_CMB_{band:03}.png"
        fig.savefig(fname_out)
        print(f"Wrote {fname_out}", flush=True)
    else:
        fname_out = f"cl_comparison_CMB_{band:03}_BB.png"
        fig.savefig(fname_out)
        print(f"Wrote {fname_out}", flush=True)
