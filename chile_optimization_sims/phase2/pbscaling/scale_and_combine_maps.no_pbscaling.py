# TODO: LAT high pass : 70 -> 30, use cosine profile,
# TODO: use r=0 for all realizations

# 2024-11-20: original version based on Reijo Keskitalo's scale_and_combine_maps.py
# 2024-11-20: revised version with fixes for beam smoothing and pixel windows of the foregrounds implemented by Shamik Ghosh
# 2024-11-27: RKe: modified to use 1-year noise simulations and separate wide/delens maps from Clem.
# 2024-12-03: RKe: modified to scale full sky noise maps to desired survey lengths

import glob
import os
import sys

import astropy.units as u
import healpy as hp
from mpi4py import MPI
import numpy as np


comm = MPI.COMM_WORLD
ntask = comm.size
rank = comm.rank

if rank == 0:
    print(f"Running with {ntask} MPI tasks")
prefix = f"{rank:04} : "

survey_lengths = [10]
survey_length_in = 1
wide_survey_max = 14  # Stop wide survey once this limit is reached
nlats = [3, 4, 5]  # Total number of LATs

# Beam FWHMs checked on 2024-11-20 against https://github.com/CMB-S4/s4mapbasedsims/blob/main/202410_s4_phase1_chile/instrument_model/cmbs4_instrument_model.tbl
fwhms = {
    "CHLAT_f020" : 9.6 * u.arcmin,
    "CHLAT_f030" : 7.8 * u.arcmin,
    "CHLAT_f040" : 5.3 * u.arcmin,
    "CHLAT_f090" : 2.1 * u.arcmin,
    "CHLAT_f150" : 1.3 * u.arcmin,
    "CHLAT_f220" : 0.95 * u.arcmin,
    "CHLAT_f280" : 0.83 * u.arcmin,
    "SAT_f030" : 79.2 * u.arcmin,
    "SAT_f040" : 56.6 * u.arcmin,
    "SAT_f085" : 23.6 * u.arcmin,
    "SAT_f095" : 21.2 * u.arcmin,
    "SAT_f145" : 15.0 * u.arcmin,
    "SAT_f155" : 13.9 * u.arcmin,
    "SAT_f220" : 9.4 * u.arcmin,
    "SAT_f280" : 7.8 * u.arcmin,
}

# Combine CMB, foregrounds and noise that is scaled to various survey lengths

noisedir_in = f"noise_1year_no_pbscaling"

fgroot = "/global/cfs/cdirs/cmbs4/chile_optimization/simulations/phase1/input_sky"
fgdir = f"{fgroot}/combined_foregrounds_mediumcomplexity_norg"

cmbdir = "/global/cfs/cdirs/cmb/data/generic/cmb/ffp10/mc/scalar"

tensordir = "/global/cfs/cdirs/cmb/data/generic/cmb/ffp10/mc/tensor"

def inv_map(m):
    """Return 1 / m"""
    minv = np.zeros_like(m)
    good = m != 0
    minv[good] = 1 / m[good]
    return minv

def sqrt_inv(m):
    """Return 1 / SQRT(m)"""
    minv = np.zeros_like(m)
    good = m != 0
    minv[good] = 1 / np.sqrt(m[good])
    return minv


def highpass(lmin, lmax):
    highpass = np.ones(lmax + 1)
    w = int(lmin / 10)
    highpass[:lmin - w] = 0
    ind = slice(lmin - w, lmin + w + 1)
    n = highpass[ind].size
    x = np.linspace(0, np.pi, n)
    highpass[ind] = (1 - np.cos(x)) / 2
    return highpass


ijob = -1
for band, fwhm in fwhms.items():
    alt_band = {
        "CHLAT_f020" : "f020",
        "CHLAT_f030" : "f026",
        "CHLAT_f040" : "f039",
        "CHLAT_f090" : "f092",
        "CHLAT_f150" : "f148",
        "CHLAT_f220" : "f227",
        "CHLAT_f280" : "f286",
        "SAT_f030" : "f026",
        "SAT_f040" : "f039",
        "SAT_f085" : "f085",
        "SAT_f095" : "f095",
        "SAT_f145" : "f145",
        "SAT_f155" : "f155",
        "SAT_f220" : "f227",
        "SAT_f280" : "f286",
    }[band]
    print(prefix + f"band = {band}, alt_band = {alt_band}, fwhm = {fwhm}")

    # Get the foregrounds

    if band.startswith("SAT"):
        fname_fg = f"{fgdir}/" \
            f"cmbs4_combined_foregrounds_mediumcomplexity_norg_uKCMB_{band}_nside512.fits"
        nside = 512
        lmax = 2 * nside
    else:
        fname_fg = f"{fgdir}/" \
            f"cmbs4_combined_foregrounds_mediumcomplexity_norg_uKCMB_{band}_nside4096.fits"
        nside = 2048
        lmax = 2 * nside

    fg = None

    for mc in range(10):

        ijob += 1
        if ijob % ntask != rank:
            continue

        # Check if the last maps for this MC were already written

        survey_length = survey_lengths[-1]
        if band.startswith("SAT"):
            fname_total = f"total_{survey_length:02}_years/" \
                f"phase2_total_{alt_band}_SAT_mc_{mc:04}.fits"
            if os.path.isfile(fname_total):
                print(prefix + f"    {fname_total} exists, skipping mc = {mc}")
                continue
            else:
                print(prefix + f"    {fname_total} does not exist.")
        else:
            fname_total5 = f"total_{survey_length:02}_years/" \
                f"phase2_total_{alt_band}_5LAT_mc_{mc:04}.fits"
            if os.path.isfile(fname_total5) and False:
                print(prefix + f"    {fname_total5} exists, skipping mc = {mc}")
                continue

        # Proceed with the simulation

        if fg is None:
            # First time foregrounds are needed
            print(prefix + f"    Reading {fname_fg}")
            fg = hp.read_map(fname_fg, None) * 1e-6  # to K_CMB
            print(prefix + f"    Expanding foreground in alm")
            fg_alms = hp.map2alm(fg, lmax=lmax, use_pixel_weights=True)

        # Process all survey lengths for this band and this MC

        # Load the CMB realization

        print(prefix + f"    MC = {mc}")
        fname_cmb = f"{cmbdir}/ffp10noaber_lensed_scl_cmb_000_alm_mc_{mc:04}.fits"
        print(prefix + f"        Reading {fname_cmb}")
        alms = hp.read_alm(fname_cmb, hdu=(1, 2, 3))
        lmax_in = hp.Alm.getlmax(alms[0].size)
        if lmax != lmax_in:
            alms = hp.resize_alm(alms, lmax_in, lmax_in, lmax, lmax)

        # Add primordial B-modes to even-numbered MC

        # if mc % 2 == 0:
        #     fname_tensor = f"{tensordir}/ffp10_ten_cmb_000_alm_mc_{mc:04}.fits"
        #     print(prefix + f"        Reading {fname_tensor}")
        #     tensor_alms = hp.read_alm(fname_tensor, hdu=(1, 2, 3))
        #     lmax_in = hp.Alm.getlmax(tensor_alms[0].size)
        #     if lmax != lmax_in:
        #         tensor_alms = hp.resize_alm(tensor_alms, lmax_in, lmax_in, lmax, lmax)
        #     r_in = 0.01  # FFP10 simulations have r=0.01
        #     r_out = 0.003  # We want r = 0.003
        #     scale = np.sqrt(r_out / r_in)  # Power spectrum ratios in map domain
        #     for i in range(3):
        #         alms[i] += scale * tensor_alms[i]
        #     r = r_out
        # else:
        fname_tensor = None
        r = 0

        # High-pass filter the CMB

        lmin = 30
        fl = np.sqrt(highpass(lmin, lmax))
        for i in range(3):
            hp.almxfl(alms[i], fl, inplace=True)

        # UPDATE 2024-11-20: Smoothing the CMB alms with the beam FWHM

        hp.smoothalm(alms, fwhm=fwhm.to_value(u.radian), inplace=True, pol=True)

        # Make a pure CMB map for reference

        print(prefix + "        Synthesizing CMB map")
        cmb = hp.alm2map(alms, nside)  # K_CMB
        if band.startswith("SAT"):
            fname_cmb = f"cmb/sat/phase2_cmb_{alt_band}_mc_{mc:04}.fits"
        else:
            fname_cmb = f"cmb/lat/phase2_cmb_{alt_band}_mc_{mc:04}.fits"
        print(prefix + f"        Writing {fname_cmb}")
        os.makedirs(os.path.dirname(fname_cmb), exist_ok=True)
        hp.write_map(
            fname_cmb,
            cmb,
            dtype=np.float32,
            coord="C",
            column_units="K_CMB",
            extra_header=[
                ("FWHM", fwhm.to_value(u.arcmin), "arcminutes"),
                ("CMB", fname_cmb),
                ("TENSOR", fname_tensor),
                ("R", r, "tensor-scalar ratio"),
                ("LMIN_CMB", lmin, "High-pass cut-off"),
            ],
            overwrite=True,
        )

        # Combine the foregrounds with the CMB and translate into a map
        # UPDATE 2024-11-20: Beam smoothing and pixel window smoothing are now not
        #                    performed on the total signal alms.
        for i in range(3):
            alms[i] += fg_alms[i]
        print(prefix + "        Synthesizing sky map")
        sky = hp.alm2map(alms, nside)  # K_CMB

        if band.startswith("SAT"):
            flavor = "sun90"
            fname_noise_in = f"{noisedir_in}/noise_sat_{flavor}_{band[-4:]}_mc_{mc:04}.fits"
            print(prefix + f"        Reading {fname_noise_in}")
            noise_in = hp.read_map(fname_noise_in, None)
            for survey_length in survey_lengths:
                scale_sat = np.sqrt(survey_length_in / survey_length)
                noise_sat = noise_in * scale_sat
                noise_sat[noise_in == 0] = hp.UNSEEN
                noisedir_out = f"noise_{survey_length:02}_years"
                os.makedirs(noisedir_out, exist_ok=True)
                # Save a copy of the pure noise map
                fname_scaled_noise_sat = f"{noisedir_out}/phase2_noise_{alt_band}_SAT_mc_{mc:04}.fits"
                args = {
                    "dtype" : np.float32,
                    "coord" : "C",
                    "column_units" : "K_CMB",
                    "extra_header" : [
                        ("NOISE", fname_noise_fullsky, "Full sky noise"),
                        ("NYEAR", survey_length),
                    ],
                    "overwrite" : True,
                }
                print(prefix + f"        Writing {fname_scaled_noise_sat}")
                hp.write_map(fname_scaled_noise_sat, noise_sat, **args)
                # Now assemble the total map
                total_sat = sky + noise_sat
                total_sat[noise_sat == hp.UNSEEN] = hp.UNSEEN
                fname_total_sat = f"total_{survey_length:02}_years/" \
                    f"phase2_total_{alt_band}_SAT_mc_{mc:04}.fits"
                os.makedirs(os.path.dirname(fname_total_sat), exist_ok=True)
                args = {
                    "dtype" : np.float32,
                    "coord" : "C",
                    "column_units" : "K_CMB",
                    "extra_header" : [
                        ("FWHM", fwhm.to_value(u.arcmin), "arcminutes"),
                        ("FOREGRND", fname_fg),
                        ("CMB", fname_cmb),
                        ("TENSOR", fname_tensor),
                        ("R", r, "tensor-scalar ratio"),
                        ("LMIN_CMB", lmin, "High-pass cut-off"),
                        ("NOISE", fname_noise_fullsky),
                    ],
                    "overwrite" : True,
                }
                print(prefix + f"        Writing {fname_total_sat}")
                hp.write_map(fname_total_sat, total_sat, **args)
        else:
            fname_noise_wide = f"{noisedir_in}/noise_lat_wide_{band[-4:]}_mc_{mc:04}.fits"
            fname_noise_tiled = f"{noisedir_in}/noise_lat_delens_tiled_{band[-4:]}_mc_{mc:04}.fits"
            fname_noise_core = f"{noisedir_in}/noise_lat_delens_{band[-4:]}_mc_{mc:04}.fits"
            print(prefix + f"        Reading {fname_noise_wide}")
            noise_wide_full = hp.read_map(fname_noise_wide, None)
            print(prefix + f"        Reading {fname_noise_tiled}")
            noise_tiled_full = hp.read_map(fname_noise_tiled, None)
            print(prefix + f"        Reading {fname_noise_core}")
            noise_core_full = hp.read_map(fname_noise_core, None)

            fname_rhits_wide = f"{noisedir_in}/rhits_lat_wide_{band[-4:]}.fits"
            fname_rhits_tiled = f"{noisedir_in}/rhits_lat_delens_tiled_{band[-4:]}.fits"
            fname_rhits_core = f"{noisedir_in}/rhits_lat_delens_{band[-4:]}.fits"

            rhits_wide_in = hp.read_map(fname_rhits_wide)
            rhits_tiled_in = hp.read_map(fname_rhits_tiled)
            rhits_core_in = hp.read_map(fname_rhits_core)

            for nlat in nlats:
                for survey_length in survey_lengths:
                    # 2 LATs dedicated to the wide survey until
                    # 14 LAT years are acquired
                    nyear_wide = min(2 * survey_length, wide_survey_max)
                    nyear_delens = nlat * survey_length - nyear_wide
                    scale_wide = np.sqrt(survey_length_in / nyear_wide)
                    scale_delens = np.sqrt(survey_length_in / nyear_delens)
                    # Scale the wide and delens maps separately and combine using
                    # inverse variance weights
                    noise_wide = noise_wide_full * scale_wide
                    noise_tiled = noise_tiled_full * scale_delens
                    rhits_wide = rhits_wide_in / scale_wide**2
                    rhits_tiled = rhits_tiled_in / scale_wide**2
                    noise_tiled = (
                        rhits_wide * noise_wide + rhits_tiled * noise_tiled
                    ) * inv_map(rhits_wide + rhits_tiled)
                    noise_tiled[noise_tiled == 0] = hp.UNSEEN
                    # Same for the core delensing strategy
                    noise_core = noise_core_full * scale_delens
                    rhits_core = rhits_core_in / scale_wide**2
                    noise_core = (
                        rhits_wide * noise_wide + rhits_core * noise_core
                    ) * inv_map(rhits_wide + rhits_core)
                    noise_core[noise_core == 0] = hp.UNSEEN
                    noisedir_out = f"noise_{survey_length:02}_years"
                    os.makedirs(noisedir_out, exist_ok=True)
                    # Save a copy of the relative hits, if theay are not
                    # already written
                    fname_rhits_tiled = f"{noisedir_out}/phase2_rhits_{alt_band}_{nlat}LAT_tiled.fits"
                    fname_rhits_core = f"{noisedir_out}/phase2_rhits_{alt_band}_{nlat}LAT_core.fits"
                    if not os.path.isfile(fname_rhits_tiled) \
                       or not os.path.isfile(fname_rhits_core):
                        args = {
                            "dtype" : np.float32,
                            "coord" : "C",
                            "column_units" : "RHITS",
                            "extra_header" : [
                                ("NWIDE", nyear_wide),
                                ("NDELENS", nyear_delens),
                                ("NYEAR", survey_length),
                            ],
                            "overwrite" : True,
                        }
                        print(prefix + f"        Writing {fname_rhits_tiled}")
                        hp.write_map(fname_rhits_tiled, rhits_tiled, **args)
                        print(prefix + f"        Writing {fname_rhits_core}")
                        hp.write_map(fname_rhits_core, rhits_core, **args)
                    # Save a copy of the pure noise map
                    fname_scaled_noise_tiled = f"{noisedir_out}/" \
                        f"phase2_noise_{alt_band}_{nlat}LAT_tiled_mc_{mc:04}.fits"
                    fname_scaled_noise_core = f"{noisedir_out}/" \
                        f"phase2_noise_{alt_band}_{nlat}LAT_core_mc_{mc:04}.fits"
                    args = {
                        "dtype" : np.float32,
                        "coord" : "C",
                        "column_units" : "K_CMB",
                        "extra_header" : [
                            ("WNOISE", fname_noise_wide, "Full sky wide noise"),
                            ("NWIDE", nyear_wide),
                            ("NDELENS", nyear_delens),
                            ("NYEAR", survey_length),
                        ],
                        "overwrite" : True,
                    }
                    args["extra_header"] += (
                        "DNOISE",
                        os.path.basename(fname_noise_tiled),
                        "Full sky delens noise"
                    )
                    print(prefix + f"        Writing {fname_scaled_noise_tiled}")
                    hp.write_map(fname_scaled_noise_tiled, noise_tiled, **args)
                    args["extra_header"][-1] = (
                        "DNOISE",
                        os.path.basename(fname_noise_core),
                        "Full sky delens noise"
                    )
                    print(prefix + f"        Writing {fname_scaled_noise_core}")
                    hp.write_map(fname_scaled_noise_core, noise_core, **args)
                    # Now assemble the total map
                    total_tiled = sky + noise_tiled
                    total_core = sky + noise_core
                    total_tiled[noise_tiled == hp.UNSEEN] = hp.UNSEEN
                    total_core[noise_core == hp.UNSEEN] = hp.UNSEEN
                    fname_total_tiled = f"total_{survey_length:02}_years/" \
                        f"phase2_total_{alt_band}_{nlat}LAT_tiled_mc_{mc:04}.fits"
                    fname_total_core = f"total_{survey_length:02}_years/" \
                        f"phase2_total_{alt_band}_{nlat}LAT_core_mc_{mc:04}.fits"
                    os.makedirs(os.path.dirname(fname_total_tiled), exist_ok=True)
                    args = {
                        "dtype" : np.float32,
                        "coord" : "C",
                        "column_units" : "K_CMB",
                        "extra_header" : [
                            ("FWHM", fwhm.to_value(u.arcmin), "arcminutes"),
                            ("FOREGRND", fname_fg),
                            ("CMB", fname_cmb),
                            ("TENSOR", fname_tensor),
                            ("R", r, "tensor-scalar ratio"),
                            ("LMIN_CMB", lmin, "High-pass cut-off"),
                            ("WNOISE", fname_noise_wide),
                        ],
                        "overwrite" : True,
                    }
                    args["extra_header"] += ("DNOISE", fname_noise_tiled, "Full sky delens noise")
                    print(prefix + f"        Writing {fname_total_tiled}")
                    hp.write_map(fname_total_tiled, total_tiled, **args)
                    args["extra_header"][-1] = ("DNOISE", fname_noise_core, "Full sky delens noise")
                    print(prefix + f"        Writing {fname_total_core}")
                    hp.write_map(fname_total_core, total_core, **args)

comm.Barrier()
if rank == 0:
    print("All done!")
comm.Barrier()