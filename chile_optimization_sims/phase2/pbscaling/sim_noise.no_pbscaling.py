# Alternative to sim_noise.sat.py and sim_noise.lat.py
# Simulates noise directly at the ab initio levels by only adding 1/ell noise

import glob
import os
import sys

import astropy.units as u
import healpy as hp
import matplotlib.pyplot as plt
from mpi4py import MPI
import numpy as np

from toast.utils import name_UID


comm = MPI.COMM_WORLD
ntask = comm.size
rank = comm.rank

if rank == 0:
    print(f"Running with {ntask} MPI tasks")
prefix = f"{rank:04} : "

params = {
    "sat" : {
        "nside" : 512,
        "ellmin" : 30,
        "flavors" : ["sun90max", "sun45_supplement"],
    },
    "lat" : {
        "nside" : 2048,
        "ellmin" : 30,
        "flavors" : ["lat_wide", "lat_delensing_max"],
    }
}

nyear = {
    "sun90max" : 10,
    "sun45_supplement" : 10,
    "lat_wide" : 14,
    "lat_delensing_max" : 16,
}

measurement_requirement = {
    # SAT TT requirements from /global/cfs/cdirs/cmbs4/chile_optimization/simulations/phase1/noise_sims/params_sat.dat
    # SAT EE/BB requirements from PBDR, Table 2-1
    # Joint band requirements added using inverse variance weights between the split bands
    "sat" : {
        # band    TTnoise TTknee TTslope EEnoise EEknee EEslope BBnoise BBknee BBslope
        "f030" : [   9.48,   159,   -4.3,    4.4,    60,   -2.0,    4.4,    60,   -1.7],
        "f040" : [  11.69,   159,   -4.3,    5.5,    60,   -2.0,    5.5,    60,   -1.7],
        "f085" : [   3.00,   159,   -4.3,    1.4,    60,   -2.0,    1.4,    60,   -1.7],
        "f090" : [   1.87,   159,   -4.3,   0.86,    60,   -2.0,   0.86,    60,   -1.7],
        "f095" : [   2.39,   159,   -4.3,    1.1,    60,   -2.0,    1.1,    60,   -1.7],
        "f145" : [   8.45,   195,   -3.7,    1.9,    65,   -2.0,    1.9,    60,   -3.0],
        "f150" : [   5.96,   195,   -3.7,    1.3,    65,   -2.0,    1.3,    60,   -3.0],
        "f155" : [   8.42,   195,   -3.7,    1.8,    65,   -2.0,    1.8,    60,   -3.0],
        "f220" : [  11.71,   195,   -3.7,    2.6,    65,   -2.0,    2.6,    60,   -3.0],
        "f280" : [  28.31,   195,   -3.7,    6.2,    65,   -2.0,    6.2,    60,   -3.0],
    },
    # LAT requirements from PBDR, Table 2-4 except for noise levels which are from Table 2-3
    # 20GHz knee/slope copied from 30GHz
    "lat" : {
        # band    TTnoise TTknee TTslope EEnoise EEknee EEslope BBnoise BBknee BBslope
        "f020" : [   11.9,   415,   -3.5,   16.7,   700,   -1.4,   16.7,   700,   -1.4],
        "f030" : [    6.5,   415,   -3.5,    9.2,   700,   -1.4,    9.2,   700,   -1.4],
        "f040" : [    3.0,   391,   -3.5,    4.2,   700,   -1.4,    4.2,   700,   -1.4],
        "f090" : [   0.45,  1932,   -3.5,   0.63,   700,   -1.4,   0.63,   700,   -1.4],
        "f150" : [   0.41,  3917,   -3.5,   0.41,   700,   -1.4,   0.41,   700,   -1.4],
        "f220" : [    1.3,  6740,   -3.5,    1.8,   700,   -1.4,    1.8,   700,   -1.4],
        "f280" : [    3.1,  6792,   -3.5,    4.3,   700,   -1.4,    4.3,   700,   -1.4],
    },
}

rootdir = "/global/cfs/cdirs/cmbs4/chile_optimization/simulations/phase2"

outdir = "noise_1year_no_pbscaling"
os.makedirs(outdir, exist_ok=True)


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


# Loop over each frequency

ijob = -1
for tele, teleparams in measurement_requirement.items():
    print(f"tele = {tele}")
    ellmin = params[tele]["ellmin"]
    nside = params[tele]["nside"]
    flavors = params[tele]["flavors"]
    npix = 12 * nside**2
    lmax = 2 * nside
    ell = np.arange(lmax + 1)
    for band, bandparams in teleparams.items():
        freq = int(band[1:])
        (
            ttnoise, ttknee, ttalpha,
            eenoise, eeknee, eealpha,
            bbnoise, bbknee, bbalpha,
        ) = bandparams
        for flavor in flavors:
            pattern = f"/global/cfs/cdirs/cmbs4/chile_optimization/simulations/" \
                + f"phase2/noise_depth/{flavor}_{band}_*years_cov.fits"
            try:
                fname_cov = sorted(glob.glob(pattern))[0]
            except:
                print(f"No matches to pattern '{pattern}'")
                sys.exit()
            for mc in range(100):
                ijob += 1
                if ijob % ntask != rank:
                    continue

                # Different seed for each frequency and MC
                rng = np.random.default_rng(name_UID(flavor + band) + mc)

                fname_full = f"{outdir}/noise_{flavor}_{band}_mc_{mc:04}.fits"

                if os.path.isfile(fname_full):
                    print(prefix + f"{fname_full} exists. Skipping.")
                    continue

                print(prefix + f"MC == {mc:04}")

                print(prefix + f"Loading {fname_cov}")
                cov = hp.read_map(fname_cov, None)
                cov = hp.ud_grade(cov, nside, power=2)
                cov *= nyear[flavor]  # Scale to one year instead of 10

                # Save a copy of the rhit map
                
                fname_rhit = f"{outdir}/rhits_{flavor}_{band}.fits"
                if not os.path.isfile(fname_rhit):
                    # Don't normalize the inverse covariance.  We need
                    # the absolute value to weight wide and delensing surveys
                    # correctly
                    invcov = invert_map(cov[0])
                    rhit = invcov # / np.amax(invcov)
                    print(prefix + f"Writing {fname_rhit}")
                    hp.write_map(
                        fname_rhit,
                        rhit, coord="C",
                        column_units="RHIT",
                        dtype=np.float32,
                        overwrite=True,
                    )

                # Start from a unit variance white noise map

                def raw_noise(cov=None):
                    m = rng.standard_normal(3 * npix).reshape(3,-1)
                    if cov is not None:
                        m *= np.sqrt(cov)
                    return m

                noisemap = raw_noise()
                print(prefix + f"Map to Alm")
                alm = hp.map2alm(raw_noise(), lmax=lmax, iter=0)

                # Add 1/ell noise

                for i, (knee, alpha) in enumerate([
                        (ttknee, ttalpha), (eeknee, eealpha), (bbknee, bbalpha)
                ]):
                    scale = np.zeros(lmax + 1)
                    # Scale the a_lm to create 1/ell noise
                    scale[1:] = np.sqrt((ell[1:] / knee)**alpha)
                    scale *= np.sqrt(highpass(ellmin, lmax))
                    hp.almxfl(alm[i], scale, inplace=True)

                print(prefix + f"Alm to Map")
                noisemap += hp.alm2map(alm, nside, lmax=lmax)
                noisemap *= np.sqrt(cov)

                # DEBUG begin
                # invcov = invert_map(cov[0])
                # rhit = invcov / np.amax(invcov)
                # rhit[rhit < 0.50] = 0
                # rhit_scale = invert_map(rhit)**.5
                # fsky = np.sum(rhit**2) / rhit.size
                # good = rhit != 0
                # cl = hp.anafast(noisemap * rhit_scale * rhit, lmax=3*nside, iter=0) / fsky
                # m0 = hp.read_map("noise_10_years/phase2_noise_f026_SAT_mc_0000.fits", None)
                # cl0 = hp.anafast(m0 * rhit_scale * rhit, lmax=3*nside, iter=0) / fsky
                # import matplotlib.pyplot as plt
                # fig = plt.figure()
                # for i in range(3):
                #     ax = fig.add_subplot(1, 3, 1 + i)
                #     ax.plot(cl0[i])
                #     ax.loglog(cl[i])
                # plt.show()
                # DEBUG end

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
