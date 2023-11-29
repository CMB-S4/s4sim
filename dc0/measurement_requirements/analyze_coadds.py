import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
# from toast.pixels_io import read_healpix

import requirements as req


#multipanel = False
multipanel = True

rootdir = "/global/cfs/cdirs/cmbs4/dc/dc0/staging/noise_sim/outputs_rk"

bands = {
    "LAT0_CHLAT" : (30, 40, 90, 150, 220, 280),
    #"LAT2_SPLAT" : (20, 30, 40, 90, 150, 220, 280),
    #"SAT1_SAT" : (95, 155, 220, 280),
    #"SAT2_SAT" : (85, 95, 145, 155, 220, 280),
    "spsat" : (30, 40, 85, 95, 145, 155, 220, 280),
    "splat" : (20, 30, 40, 90, 150, 220, 280),
    #"SAT3_SAT" : (30, 40, 85, 145),
}

alt_bands = {
    "chlat" : {30 : 27, 40 : 39, 90 : 93, 150 : 145, 220 : 225, 280 : 278},
    "splat" : {20 : 20, 30 : 27, 40 : 39, 90 : 93, 150 : 145, 220 : 225, 280 : 278},
    "spsat" : {30 : 30, 40 : 40, 85 : 85, 95 : 95, 145 : 145, 155 : 155, 220 : 220, 280 : 270},
}

fskies = {
    "chlat" : 0.60,
    "splat" : 0.03,
    "spsat" : 0.03,
}

for TELE in bands:
    tele = TELE.split("_")[-1].lower()
    if tele == "sat":
        tele = "spsat"
    if tele == "spsat":
        nside = 512
    else:
        nside = 4096
    npix = 12 * nside**2
    fsky_req = fskies[tele]
    lmax = 2 * nside
    ell = np.arange(lmax + 1)
    for band in bands[TELE]:
        fname_map = f"{rootdir}/coadd/{TELE}/coadd_{TELE}_f{band:03}_001of001_map.fits"
        if not os.path.isfile(fname_map):
            print(f"Not found: {fname_map}")
            continue
        fname_mask = f"/global/cfs/cdirs/cmbs4/dc/dc0/masks/mask_mediumcomplexity.{tele}.f{band:03}.fits"
        fname_cov = f"{rootdir}/coadd/{TELE}/coadd_{TELE}_f{band:03}_001of001_cov.fits"
        fname_cl = f"outputs/cl/{TELE}/coadd_{TELE}_f{band:03}_001of001_cl.fits"

        freq = alt_bands[tele][band]
        req_ell = req.ells[:lmax + 1]
        if tele == "chlat":
            fwhm = req.Chile_LAT[freq][0]
            nltt = req.NlTT_Chile_LAT[freq][:lmax + 1]
            nlee = req.NlEE_Chile_LAT[freq][:lmax + 1]
        elif tele == "splat":
            fwhm = req.Pole_LAT[freq][0]
            nltt = req.NlTT_Pole_LAT[freq][:lmax + 1]
            nlee = req.NlEE_Pole_LAT[freq][:lmax + 1]
        else:
            fwhm = req.Pole_SAT[freq][0]
            nltt = req.NlTT_Pole_SAT[freq][:lmax + 1]
            nlee = req.NlEE_Pole_SAT[freq][:lmax + 1]
        bl = req.get_bl(fwhm, ell)
        bl = 1 / bl[:lmax + 1]

        fname_cl_cmb_in = f"../multimap_sim/outputs/cl/{TELE}/input_cmb_{TELE}_f{band:03}_cl.fits"
        fname_cl_cmb_out = f"../multimap_sim/outputs/cl/{TELE}/coadd_cmb_{TELE}_f{band:03}_cl.fits"
        print(f"Loading {fname_cl_cmb_in}", flush=True)
        cl_in = hp.read_cl(fname_cl_cmb_in)
        print(f"Loading {fname_cl_cmb_out}", flush=True)
        cl_out = hp.read_cl(fname_cl_cmb_out)
        tf = cl_out[:3] / cl_in[:3]
        tf[:, :2] = 1
        # Overwrite transfer function where there is no signal to measure
        tf = tf * bl + np.ones_like(tf) * (1 - bl)
        # Save the TF
        # fname_cl = f"outputs/cl/{TELE}/tf_{TELE}_f{band:03}_001of001_cl.fits"
        fname_tf = f"outputs/cl/tf_{TELE}_f{band:03}_001of001_cl.fits"
        hp.write_cl(fname_tf, tf, overwrite=True)
        print(f"Wrote transfer function to {fname_tf}")

        if os.path.isfile(fname_cl):
            # if False:
            print(f"Loading {fname_cl}", flush=True)
            cl = hp.read_cl(fname_cl)
        else:
            outroot = os.path.dirname(fname_cl)
            os.makedirs(outroot, exist_ok=True)
            print(f"Loading {fname_map}")
            m = hp.read_map(fname_map, None)

            # Read and invert the processing mask
            print(f"Loading {fname_mask}")
            mask = np.logical_not(hp.read_map(fname_mask, dtype=bool))

            # Discard 1% of the noisiest pixels
            print(f"Loading {fname_cov}")
            w = hp.read_map(fname_cov)
            good = w > 0
            sorted_w = np.sort(w[good])
            # ngood = sorted_w.size
            # lim = sorted_w[int(.99 * ngood)]
            lim = sorted_w[int(npix * fsky_req)]
            mask[w == 0] = False
            mask[w > lim] = False

            print("Measuring C_ell")
            m[0] = hp.remove_dipole(m[0] * mask, bad=0)
            cl = hp.anafast(m * mask, lmax=lmax, iter=0)
            fsky = np.sum(mask) / mask.size
            cl /= fsky
            print(f"Writing {fname_cl}", flush=True)
            hp.write_cl(fname_cl, cl, overwrite=True)

        if multipanel:
            nrow, ncol = 1, 3
        else:
            nrow, ncol = 1, 1
        fig = plt.figure(figsize=[4 * ncol, 3 * nrow])
        ell = np.arange(lmax + 1)
        ellnorm = ell * (ell + 1) / (2 * np.pi) * 1e12

        iplot = 1
        ax = fig.add_subplot(nrow, ncol, iplot)
        ax.set_title(f"TT {band}GHz")
        ax.set_xlabel("Multipole, $\ell$")
        ax.set_ylabel("D$\ell$ [$\mu$K$^2$]")
        fiducial_ell = req.fiducial_ell[:lmax + 1]
        fiducial_TT = req.fiducial_TT[:lmax + 1]
        ax.loglog(fiducial_ell, fiducial_TT, "k", label="CMB")
        mr = {
            "LAT0_CHLAT" : "MR 2.0 & 3.1",
            "spsat" : None,
            "splat" : "MR 3.2",
        }[TELE]
        ax.loglog(req_ell, nltt, label=mr)
        ax.loglog(ell[2:], (ellnorm * cl[0])[2:], label=f"DC0")
        ax.loglog(ell[2:], (ellnorm * cl[0] / bl / tf[0])[2:], label=f"TF-corrected DC0")
        ax.legend(loc="best")
        ax.set_ylim([np.amin(fiducial_TT) * 1e-5, np.amax(fiducial_TT) * 1e5])
        if not multipanel:
            fname_out = f"cl_comparison_noise_{tele}_{band:03}_TT.png"
            fig.savefig(fname_out)
            print(f"Wrote {fname_out}", flush=True)
            fname_out = f"cl_comparison_noise_{tele}_{band:03}_TT.pdf"
            fig.savefig(fname_out)
            print(f"Wrote {fname_out}", flush=True)
            plt.clf()
            iplot = 0

        iplot += 1
        ax = fig.add_subplot(nrow, ncol, iplot)
        ax.set_title(f"EE {band}GHz")
        ax.set_xlabel("Multipole, $\ell$")
        ax.set_ylabel("D$\ell$ [$\mu$K$^2$]")
        fiducial_EE = req.fiducial_EE[:lmax + 1]
        ax.loglog(fiducial_ell, fiducial_EE, "k", label="CMB")
        mr = {
            "LAT0_CHLAT" : "MR2 2.0",
            "spsat" : "MR 1.1",
            "splat" : "MR 1.2",
        }[TELE]
        ax.loglog(req_ell, nlee, label=mr)
        ax.loglog(ell[2:], (ellnorm * cl[1])[2:], label=f"DC0")
        ax.loglog(ell[2:], (ellnorm * cl[1] / bl / tf[1])[2:], label=f"TF-corrected DC0")
        ax.legend(loc="best")
        ax.set_ylim([np.amin(fiducial_EE) * 1e-5, np.amax(fiducial_EE) * 1e5])
        if not multipanel:
            fname_out = f"cl_comparison_noise_{tele}_{band:03}_EE.png"
            fig.savefig(fname_out)
            print(f"Wrote {fname_out}", flush=True)
            fname_out = f"cl_comparison_noise_{tele}_{band:03}_EE.pdf"
            fig.savefig(fname_out)
            print(f"Wrote {fname_out}", flush=True)
            plt.clf()
            iplot = 0

        iplot += 1
        ax = fig.add_subplot(nrow, ncol, iplot)
        ax.set_title(f"BB {band}GHz")
        ax.set_xlabel("Multipole, $\ell$")
        ax.set_ylabel("D$\ell$ [$\mu$K$^2$]")
        fiducial_BB = req.fiducial_BB[:lmax + 1]
        ax.loglog(fiducial_ell, fiducial_BB, "k", label="CMB")
        ax.loglog(req_ell, nlee, label="MR 2.0")
        ax.loglog(ell[2:], (ellnorm * cl[2])[2:], label=f"DC0")
        ax.loglog(ell[2:], (ellnorm * cl[2] / bl / tf[2])[2:], label=f"DC0")
        ax.set_ylim([np.amin(fiducial_BB) * 1e-5, np.amax(fiducial_BB) * 1e5])
        fig.tight_layout()
        if multipanel: 
            fname_out = f"cl_comparison_noise_{tele}_{band:03}.png"
            fig.savefig(fname_out)
            print(f"Wrote {fname_out}", flush=True)
            fname_out = f"cl_comparison_noise_{tele}_{band:03}.pdf"
            fig.savefig(fname_out)
            print(f"Wrote {fname_out}", flush=True)
        else:
            fname_out = f"cl_comparison_noise_{tele}_{band:03}_BB.png"
            fig.savefig(fname_out)
            print(f"Wrote {fname_out}", flush=True)
            fname_out = f"cl_comparison_noise_{tele}_{band:03}_BB.pdf"
            fig.savefig(fname_out)
            print(f"Wrote {fname_out}", flush=True)
