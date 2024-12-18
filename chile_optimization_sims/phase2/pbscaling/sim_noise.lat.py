# Based on https://bicep.rc.fas.harvard.edu/CMB-S4/analysis_logbook/20241202_dc12_lat_1yr

import os
import sys

import astropy.units as u
import healpy as hp
import matplotlib.pyplot as plt
from mpi4py import MPI
import numpy as np


comm = MPI.COMM_WORLD
ntask = comm.size
rank = comm.rank

if rank == 0:
    print(f"Running with {ntask} MPI tasks")
prefix = f"{rank:04} : "

rootdir1 = "/global/cfs/cdirs/cmbs4/chile_optimization/simulations/phase1"
rootdir2 = "/global/cfs/cdirs/cmbs4/chile_optimization/simulations/phase2"

outdir_fullsky = "noise_fullsky"
os.makedirs(outdir_fullsky, exist_ok=True)

outdir_rhits = "rhits"
os.makedirs(outdir_rhits, exist_ok=True)


bc_to_band = {
    20.00 : "f020",
    25.75 : "f030",
    38.75 : "f040",
    91.50 : "f090",
    148.5 : "f150",
    227.0 : "f220",
    285.5 : "f280",
}

arr = np.genfromtxt(
    "noise_sims/params_lat.dat",
    comments="%",
    names=[
        "BC", "BW", "FWHM",
        "TTnoise", "TTknee", "TTalpha",
        "EEnoise", "EEknee", "EEalpha",
        "BBnoise", "BBknee", "BBalpha",
        "ellmin", "nside", "rseed",
    ],
)

reldetyrs_delens = np.genfromtxt(
    "noise_sims/reldetyrs_phase1_delens_1chlatyr.dat",
    names=["band", "years"],
)
reldetyrs_wide = np.genfromtxt(
    "noise_sims/reldetyrs_phase1_wide_1chlatyr.dat",
    names=["band", "years"],
)

def invert_map(m):
    """Invert non-zero pixels"""
    result = np.zeros_like(m)
    good = m != 0
    result[good] = 1 / m[good]
    return result

# Load the SPLAT relative hit map, used as the basis of the LAT noise scaling

rhit0 = hp.read_map("noise_sims/alt1_hits_splat_n2048.fits", None)
fsky = np.sum(rhit0**2) / rhit0.size
good = rhit0 != 0
rhit0_scale = invert_map(rhit0)**.5

# ref_map = hp.read_map("/global/cfs/cdirs/cmbs4/awg/lowellbb/expt_xx/12a1lat/noise/map/noise_f020_b11.40_ellmin70_map_2048_mc_0000.fits", None)
#ref_cl = hp.anafast(ref_map * rhit0_scale * rhit0, lmax=4096, iter=0) / fsky
# ref_cl = hp.anafast(ref_map, lmax=4096, iter=0)

# Loop over each frequency

ijob = -1
for irow, row in enumerate(arr):
    bc, bw, fwhm, ttnoise, ttknee, ttalpha, eenoise, eeknee, eealpha, bbnoise, bbknee, bbalpha, ellmin, nside, rseed = row
    band = bc_to_band[bc]
    ellmin = int(ellmin)
    nside = int(nside)
    npix = 12 * nside**2
    lmax = 2 * nside
    ell = np.arange(lmax + 1)

    fname_rhit_wide = f"{outdir_rhits}/rhits_wide_{band}.fits"
    fname_rhit_delens_tiled = f"{outdir_rhits}/rhits_delens_tiled_{band}.fits"
    fname_rhit_delens_core = f"{outdir_rhits}/rhits_delens_core_{band}.fits"
    if rank == 0 and (
            not os.path.isfile(fname_rhit_wide) \
            or not os.path.isfile(fname_rhit_delens_tiled) \
            or not os.path.isfile(fname_rhit_delens_core)
    ):

        # We'll use Colin's reldetyrs from Phase1 so we need the relative
        # survey weight between Phase1 and Phase2

        cov0_wide = hp.read_map(f"{rootdir1}/noise_depth/lat_wide_{band}_14years_cov.fits")
        cov_wide = hp.read_map(f"{rootdir2}/noise_depth/lat_wide_{band}_14years_cov.fits")
        invcov0_wide = invert_map(cov0_wide)
        invcov_wide = invert_map(cov_wide)
        relative_weight_wide = np.sum(invcov_wide) / np.sum(invcov0_wide)

        cov0_delens = hp.read_map(f"{rootdir1}/noise_depth/lat_delensing_{band}_16years_cov.fits")
        invcov0_delens = invert_map(cov0_delens)

        cov_delens_tiled = hp.read_map(f"{rootdir2}/noise_depth/lat_delensing_tiled_{band}_16years_cov.fits")
        invcov_delens_tiled = invert_map(cov_delens_tiled)
        relative_weight_delens_tiled = np.sum(invcov_delens_tiled) / np.sum(invcov0_delens)

        cov_delens_core = hp.read_map(f"{rootdir2}/noise_depth/lat_delensing_core_{band}_16years_cov.fits")
        invcov_delens_core = invert_map(cov_delens_core)
        relative_weight_delens_core = np.sum(invcov_delens_core) / np.sum(invcov0_delens)

        # Colin's relative detector years from Phase-1 scale the results to
        # a single observing year of a single LAT

        rel_years_wide = reldetyrs_wide[irow]["years"]
        rel_years_delens = reldetyrs_delens[irow]["years"]

        # Scale from SPLAT to Phase-1 and to Phase-2

        def invcov_to_rhit(invcov):
            invcov = hp.ud_grade(invcov, 2048)
            return invcov * np.sum(rhit0) / np.sum(invcov)

        rhit_wide = invcov_to_rhit(invcov_wide) \
            * rel_years_wide * relative_weight_wide
        rhit_delens_tiled = invcov_to_rhit(invcov_delens_tiled) \
            * rel_years_delens * relative_weight_delens_tiled
        rhit_delens_core = invcov_to_rhit(invcov_delens_core) \
            * rel_years_delens * relative_weight_delens_core

        # Save the relative hitmaps for future

        args = {
            "coord" : "C",
            "column_units" : "RHIT",
            "dtype" : np.float32,
            "overwrite" : True,
            "extra_header" : [
                ("nyear", 1, "Number of observing years"),
                None,
            ],
        }
        args["extra_header"][-1] = (
            "rweight", "relative_weight_wide", "Phase2/Phase1 survey weight"
        )
        hp.write_map(fname_rhit_wide, rhit_wide, **args)
        print(prefix + f"Wrote {fname_rhit_wide}")
        args["extra_header"][-1] = (
            "rweight", "relative_weight_delens_tiled", "Phase2/Phase1 survey weight"
        )
        hp.write_map(fname_rhit_delens_tiled, rhit_delens_tiled, **args)
        print(prefix + f"Wrote {fname_rhit_delens_tiled}")
        args["extra_header"][-1] = (
            "rweight", "relative_weight_delens_core", "Phase2/Phase1 survey weight"
        )
        hp.write_map(fname_rhit_delens_core, rhit_delens_core, **args)
        print(prefix + f"Wrote {fname_rhit_delens_core}")

    for mc in range(10):
        ijob += 1
        if ijob % ntask != rank:
            continue

        # Different seed for each frequency and MC
        np.random.seed(347164 + int(bc) * 1000 + mc)
    
        fname_full = f"{outdir_fullsky}/noise_lat_full_sky_{band}_mc_{mc:04}.fits"

        if os.path.isfile(fname_full):
            print(prefix + f"{fname_full} exists. Skipping.")
            continue

        print(prefix + f"MC == {mc:04}")

        # Start from a unit variance white noise map

        noisemap_in = np.random.randn(3 * npix).reshape(3,-1)
        print(prefix + f"Map to Alm")
        alm = hp.map2alm(noisemap_in, lmax=lmax, iter=0)

        # Measure weighted noise in the patch to determine appropriate scaling

        print(prefix + f"Anafast")
        cl = hp.anafast(noisemap_in * rhit0_scale * rhit0, lmax=lmax, iter=0) / fsky
        for i, noise in enumerate([ttnoise, eenoise, bbnoise]):
            noiselevel = (noise * 1e-6 / 60 * np.pi / 180)**2  # uK.arcmin -> (K.rad)**2
            # This scaling will recover the prescribed noise level
            scale = np.sqrt(noiselevel / np.mean(cl[i, 2:]))
            alm[i] *= scale

        # Add 1/ell noise and high-pass

        for i, (knee, alpha) in enumerate([
                (ttknee, ttalpha), (eeknee, eealpha), (bbknee, bbalpha)
        ]):
            scale = np.zeros(lmax + 1)
            scale[ellmin:] = np.sqrt(1 + (ell[ellmin:] / knee)**alpha)
            hp.almxfl(alm[i], scale, inplace=True)

        print(prefix + f"Alm to Map")
        noisemap_full = hp.alm2map(alm, nside, lmax=lmax)
        # cl_out = hp.anafast(noisemap_full * rhit0_scale * rhit0, lmax=lmax, iter=0) / fsky
        # cl_out = hp.anafast(noisemap_full, lmax=lmax, iter=0)

        hp.write_map(
            fname_full,
            noisemap_full,
            coord="C",
            column_units="K_CMB",
            dtype=np.float32,
            overwrite=True,
        )
        print(prefix + f"Wrote {fname_full}")
