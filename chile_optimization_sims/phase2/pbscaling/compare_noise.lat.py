import glob
import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


recompute = False


def invert_map(m):
    """Invert non-zero pixels"""
    result = np.zeros_like(m)
    good = m != 0
    result[good] = 1 / m[good]
    return result


def get_cl(fname_noise, fname_rhit, lmax, recompute):
    fname_cl = fname_noise.replace(".fits", "_cl.fits")
    if fname_noise == fname_cl:
        raise RuntimeError("Failed to synthesize C_ell filename")
    if os.path.isfile(fname_cl) and not recompute:
        print(f"Loading {fname_cl}")
        cl = hp.read_cl(fname_cl)
    else:
        print(f"Loading {fname_rhit}")
        rhit = hp.read_map(fname_rhit)
        rhit /= np.amax(rhit)  # ensure peak normalization
        # Use best 3% for power spectrum
        limit = np.sort(rhit)[int(rhit.size * .97)]
        rhit[rhit < limit] = 0  # noise weights
        # rhit = rhit > limit  # binary mask
        fsky = np.sum(rhit**2) / rhit.size
        print(f"Loading {fname_noise}")
        noise = hp.read_map(fname_noise, None)
        print(f"Running Anafast on {fname_noise}")
        cl = hp.anafast(rhit * noise, lmax=lmax, iter=0) / fsky
        print(f"Writing {fname_cl}")
        hp.write_cl(fname_cl, cl, overwrite=True)
    if np.any(cl > 1):
        msg = f"Anomalous cl in {fname_cl}. Check {fname_noise}"
        raise RuntimeError(msg)
    return cl


def noise_model(ell, level, knee, alpha):
    return level * (1 + (ell / knee)**alpha)


def fit_noise(cl):
    # fit a noise model
    ell = np.arange(cl.size)
    ind = slice(100, int(2 * cl.size / 3))
    p0 = (np.mean(cl[ind]), 700, -1)
    popt, pcov = curve_fit(noise_model, ell[ind], cl[ind], p0=p0)
    return popt


nside = 2048
lmax = 3 * nside
ell = np.arange(lmax + 1)

bands = {
    "f020" : "f020",
    "f026" : "f030",
    "f039" : "f040",
    "f092" : "f090",
    "f148" : "f150",
    "f227" : "f220",
    "f286" : "f280",
}

for nlat in 3, 4, 5:

    nrow, ncol = 2, 4
    fig = plt.figure(figsize=[ncol * 4, nrow * 4])
    fig.suptitle(f"Delensing survey, {nlat}LAT")
    iplot = 0
    for alt_band, band in bands.items():
        iplot += 1
        print(f"{alt_band} {band}")
        fname_rhit = f"rhits/rhits_delens_{band}.fits"
        fname_rhit0 = "/global/cfs/cdirs/cmbs4/chile_optimization/simulations/" \
            f"phase1/noise_sims/phase1_hits_{nlat}chlat.fits"

        nyear = 10
        fname0 = f"/global/cfs/cdirs/cmbs4/chile_optimization/simulations/" \
            f"phase1/noise_{nyear:02}_years/phase1_noise_{alt_band}_{nlat}LAT_mc_0000.fits"
        fname1 = f"with_pbscaling/noise_{nyear:02}_years/" \
            f"phase2_noise_{alt_band}_{nlat}LAT_mc_0000.fits"
        fname2 = f"no_pbscaling_no_artifact/noise_{nyear:02}_years/" \
            f"phase2_noise_{alt_band}_{nlat}LAT_mc_0000.fits"

        cl0 = get_cl(fname0, fname_rhit0, lmax, recompute)
        cl1 = get_cl(fname1, fname_rhit, lmax, recompute)
        cl2 = get_cl(fname2, fname_rhit, lmax, recompute)

        params0 = fit_noise(cl0[2])
        params1 = fit_noise(cl1[2])
        params2 = fit_noise(cl2[2])

        level0 = params0[0]
        level1 = params1[0]
        level2 = params2[0]

        ratio10 = level1 / level0
        ratio21 = level2 / level1

        depth0 = np.sqrt(level0) * 1e6 * 180 / np.pi * 60
        depth1 = np.sqrt(level1) * 1e6 * 180 / np.pi * 60
        depth2 = np.sqrt(level2) * 1e6 * 180 / np.pi * 60

        ax = fig.add_subplot(nrow, ncol, iplot)
        ax.set_title(
            f"{band} : " + r"C$_\ell^\mathrm{BB}$"
            + f" ratios = {ratio10:.3f}"
            + f", {ratio21:.3f}"
        )
        ind = slice(2, lmax + 1)
        ax.loglog(ell[ind], cl0[2][ind], label=f"Phase 1 : depth = {depth0:.3f}", color="tab:blue")
        ax.loglog(ell[ind], noise_model(ell[ind], *params0), "--", color="tab:blue")
        ax.loglog(ell[ind], cl1[2][ind], label=f"Phase 2 : depth = {depth1:.3f}", color="tab:orange")
        ax.loglog(ell[ind], noise_model(ell[ind], *params1), "--", color="tab:orange")
        ax.loglog(ell[ind], cl2[2][ind], label=f"Phase 2 : depth = {depth2:.3f} (no BK)", color="tab:green")
        ax.loglog(ell[ind], noise_model(ell[ind], *params2), "--", color="tab:green")
        ax.legend(loc="best")

    fig.tight_layout()
    fig.savefig(f"lat_noise_comparison.{nlat}LAT.png")
