import glob
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
    25.00 : "f030",
    40.00 : "f040",
    85.00 : "f085",
    90.00 : "f090",
    95.00 : "f095",
    145.0 : "f145",
    150.0 : "f150",
    155.0 : "f155",
    220.0 : "f220",
    280.0 : "f280",
}

arr = np.genfromtxt(
    rootdir2 + "/noise_sims/params_sat.dat",
    comments="%",
    names=[
        "BC", "BW", "FWHM",
        "TTnoise", "TTknee", "TTalpha",
        "EEnoise", "EEknee", "EEalpha",
        "BBnoise", "BBknee", "BBalpha",
        "ellmin", "nside",
    ],
)

reldetyrs = np.genfromtxt(
    rootdir2 + "/noise_sims/reldetyrs_phase1_chsat.dat",
    comments="%",
    names=["band", "years"],
)


def invert_map(m):
    """Invert non-zero pixels"""
    result = np.zeros_like(m)
    good = m != 0
    result[good] = 1 / m[good]
    return result


def highpass(lmin, lmax):
    highpass = np.ones(lmax + 1)
    w = int(lmin / 10)
    highpass[:lmin - w] = 0
    ind = slice(lmin - w, lmin + w + 1)
    n = highpass[ind].size
    x = np.linspace(0, np.pi, n)
    highpass[ind] = (1 - np.cos(x)) / 2
    return highpass


# Load the SPSAT relative hit map, used as the basis of the SAT noise scaling

rhit0 = hp.read_map(rootdir2 + "/noise_sims/alt1_hits_sat.fits", None)
fsky = np.sum(rhit0**2) / rhit0.size
good = rhit0 != 0
rhit0_scale = invert_map(rhit0)**.5

# ref_map = hp.read_map("/global/cfs/cdirs/cmbs4/awg/lowellbb/expt_xx/12a1lat/noise/map/noise_f020_b11.40_ellmin70_map_2048_mc_0000.fits", None)
#ref_cl = hp.anafast(ref_map * rhit0_scale * rhit0, lmax=4096, iter=0) / fsky
# ref_cl = hp.anafast(ref_map, lmax=4096, iter=0)

# Loop over each frequency

ijob = -1
for irow, row in enumerate(arr):
    (
        bc, bw, fwhm,
        ttnoise, ttknee, ttalpha,
        eenoise, eeknee, eealpha,
        bbnoise, bbknee, bbalpha,
        ellmin, nside,
    ) = row
    band = bc_to_band[bc]
    ellmin = int(ellmin)
    nside = int(nside)
    npix = 12 * nside**2
    lmax = 2 * nside
    ell = np.arange(lmax + 1)

    # Joint bands do not have a phase-1 reference so just use scale from the
    # nearest available band
    if band == "f090":
        ref_band = "f085"
    elif band == "f150":
        ref_band = "f145"
    else:
        ref_band = band

    fname_rhit = f"{outdir_rhits}/rhits_sat_{band}.fits"
    if rank == 0 and not os.path.isfile(fname_rhit):
        # We'll use Colin's reldetyrs from Phase1 so we need the relative
        # survey weight between Phase1 and Phase2

        fname0 = glob.glob(f"{rootdir1}/noise_depth/sat_{ref_band}_*years_cov.fits")[0]
        cov0 = hp.read_map(fname0)
        fname = glob.glob(f"{rootdir2}/noise_depth/sun90max_{band}_*years_cov.fits")[0]
        cov = hp.read_map(fname)
        invcov0 = invert_map(cov0)
        invcov = invert_map(cov)
        relative_weight = np.sum(invcov) / np.sum(invcov0)

        # Colin's relative detector years from Phase-1 scale the results to
        # 10 years of 9 SATs.

        rel_years = reldetyrs[irow]["years"]

        # Scale from SPSAT to Phase-1 and to Phase-2

        def invcov_to_rhit(invcov):
            # invcov = hp.ud_grade(invcov, 512)
            return invcov * np.sum(rhit0) / np.sum(invcov)

        rhit = invcov_to_rhit(invcov) * rel_years * relative_weight

        # Scale from 10 years to one year

        rhit /= 10

        # Save the relative hitmaps for future

        args = {
            "coord" : "C",
            "column_units" : "RHIT",
            "dtype" : np.float32,
            "overwrite" : True,
            "extra_header" : [
                ("nyear", 1, "Number of observing years"),
                ("rweight", relative_weight, "Phase2/Phase1 survey weight"),
            ]
        }
        hp.write_map(fname_rhit, rhit, **args)
        print(prefix + f"Wrote {fname_rhit}")

    for mc in range(3):
        ijob += 1
        if ijob % ntask != rank:
            continue

        # Different seed for each frequency and MC
        rng = np.random.default_rng(93654635812 + int(bc) * 10000 + mc)
    
        fname_full = f"{outdir_fullsky}/noise_sat_full_sky_{band}_mc_{mc:04}.fits"

        if os.path.isfile(fname_full):
            print(prefix + f"{fname_full} exists. Skipping.")
            continue

        print(prefix + f"MC == {mc:04}")

        # Start from a unit variance white noise map

        def raw_noise():
            return rng.standard_normal(3 * npix).reshape(3,-1)

        noisemap = raw_noise()
        print(prefix + f"Map to Alm")
        alm = hp.map2alm(raw_noise(), lmax=lmax, iter=0)

        # Measure weighted noise in the patch to determine appropriate scaling

        print(prefix + f"Anafast")
        cl = hp.anafast(noisemap * rhit0_scale * rhit0, lmax=lmax, iter=0) / fsky
        scales = []
        for i, noise in enumerate([ttnoise, eenoise, bbnoise]):
            noiselevel = (noise * 1e-6 / 60 * np.pi / 180)**2  # uK.arcmin -> (K.rad)**2
            # This scaling will recover the prescribed noise level
            scale = np.sqrt(noiselevel / np.mean(cl[i, 2:]))
            alm[i] *= scale
            scales.append(scale)
        noisemap[0] *= scales[0]
        noisemap[1:] *= scales[2]

        # Add 1/ell noise and high-pass

        for i, (knee, alpha) in enumerate([
                (ttknee, ttalpha), (eeknee, eealpha), (bbknee, bbalpha)
        ]):
            # Scale the a_lm to create 1/ell noise
            scale = np.zeros(lmax + 1)
            scale[1:] = np.sqrt((ell[1:] / knee)**alpha)
            scale *= np.sqrt(highpass(ellmin, lmax))
            hp.almxfl(alm[i], scale, inplace=True)

        print(prefix + f"Alm to Map")
        noisemap += hp.alm2map(alm, nside, lmax=lmax)

        hp.write_map(
            fname_full,
            noisemap,
            coord="C",
            column_units="K_CMB",
            dtype=np.float32,
            overwrite=True,
        )
        print(prefix + f"Wrote {fname_full}")

comm.Barrier()
if rank == 0:
    print("All done!")
comm.Barrier()
