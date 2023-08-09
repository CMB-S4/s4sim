import os
import sys

import numpy as np
import healpy as hp
from toast.pixels_io_healpix import write_healpix


complexity = "medium"
inroot = "/pscratch/sd/z/zonca/cmbs4/202305_dc0"

outdir = "input_maps"
os.makedirs(outdir, exist_ok=True)

bands = {
    "CHLAT" : ("chlat", 4096, [30, 40, 90, 150, 220, 280]),
    "SPLAT" : ("splat", 4096, [20, 30, 40, 90, 150, 220, 280]),
    "SAT" : ("spsat", 512, [30, 40, 85, 95, 145, 155, 220, 280]),
}

snr_limits_fg = {"CHLAT" : 1000, "SPLAT" : 20000, "SAT" : 50000}
snr_limits_radio = {"CHLAT" : 30, "SPLAT" : 1000, "SAT" : 10000}
mask_fwhm = np.radians(1)
mask_lmax = 512
gal_mask_frac = 0.1

measurement_requirements = {
    #  measurement requirement in uK.arcmin, PBDR, Table 2-1, 2-2, 2-3
    "CHLAT" : {  # T / P noise
        30 :  (21.8, 30.23),
        40 :  (12.4, 16.53),
        90 :  ( 2.0,  2.68),
        150 : ( 2.0,  2.96),
        220 : ( 6.9,  9.78),
        280 : (16.7, 23.93),
    },
    "SPLAT" : {  # T / P noise
        20 :  (9.4, 13.16),
        30 :  (4.6,  6.5),
        40 :  (3.0,  4.15),
        90 :  (0.45, 0.63),
        150 : (0.41, 0.59),
        220 : (1.3,  1.83),
        280 : (3.1,  4.34),
    },
    "SAT" : {  # E / B noise
        30 :  (3.7,  3.5),
        40 :  (4.7,  4.5),
        85 :  (0.93, 0.88),
        95 :  (0.82, 0.78),
        145 : (1.3,  1.2),
        155 : (1.33, 1.3),
        220 : (3.5,  3.5),
        280 : (8.1,  6.0),
    },
}

for TELESCOPE, (telescope, nside, freqs) in bands.items():
    # Start by assembling a joint Galactic mask for this telescope
    gal_mask = None
    for freq in freqs:
        fname_in = os.path.join(
            inroot,
            f"combined_foregrounds_{complexity}complexity",
            f"cmbs4_combined_foregrounds_mediumcomplexity_uKCMB_{TELESCOPE}_f{freq:03}_nside{nside}.fits",
        )
        if not os.path.isfile(fname_in):
            raise RuntimeError(f"Input file does not exist: {fname_in}")
        print(f"Reading {fname_in}")
        m = hp.read_map(
            fname_in,
            0,
            nest=True,
            verbose=False,
            dtype=np.float32,
        )
        limit = np.percentile(m, 100 * (1 - gal_mask_frac))
        freq_mask = m > limit
        fname_out = f"{outdir}/galmask_{freq:03}_{int(gal_mask_frac * 100):02}pc.fits"
        write_healpix(fname_out, freq_mask, coord="C", nest=True, overwrite=True)
        print(f"Wrote {fname_out}")
        if gal_mask is None:
            gal_mask = freq_mask
        else:
            gal_mask[freq_mask] = True

    gal_mask = hp.smoothing(gal_mask, lmax=mask_lmax, fwhm=mask_fwhm, nest=True) > 0.25
    fname_out = f"{outdir}/galmask_{int(gal_mask_frac * 100):02}pc.fits"
    write_healpix(fname_out, gal_mask, coord="C", nest=True, overwrite=True)
    print(f"Wrote {fname_out}")

    for freq in freqs:
        fname_out = os.path.join(outdir, f"foreground_mediumcomplexity.{telescope}.f{freq:03}.h5")
        if os.path.isfile(fname_out) and False:
            print(f"Output file exists: {fname_out}")
            continue
        fname_in = os.path.join(
            inroot,
            f"combined_foregrounds_{complexity}complexity",
            f"cmbs4_combined_foregrounds_mediumcomplexity_uKCMB_{TELESCOPE}_f{freq:03}_nside{nside}.fits",
        )
        if not os.path.isfile(fname_in):
            raise RuntimeError(f"Input file does not exist: {fname_in}")
        print(f"Reading {fname_in}")
        m = hp.read_map(
            fname_in,
            None,
            nest=True,
            verbose=False,
            dtype=np.float32,
        )
        write_healpix(fname_out, m * 1e-6, coord="C", nest=True, overwrite=True)
        print(f"Wrote {fname_out}")

        # Build processing masks
        pixarea = hp.nside2pixarea(nside, degrees=True)
        mr_T, mr_P = measurement_requirements[TELESCOPE][freq]
        if TELESCOPE == "SAT":
            requirement = np.sqrt(2) * max(mr_T, mr_P)
        else:
            requirement = mr_T
        sdev = requirement / 60 / np.sqrt(pixarea)  # Translate to pixel RMS
        mask = gal_mask.copy()
        print(f"  Mask fraction after fg : {np.sum(mask)/mask.size:.3f}")
        # mask radio sources separately
        fname_radio = os.path.join(
            inroot,
            f"radio_rg1",
            f"cmbs4_radio_rg1_uKCMB_{TELESCOPE}_f{freq:03}_nside{nside}.fits",
        )
        radio = hp.read_map(fname_radio, nest=True)
        snr_radio = (radio / sdev)**2
        mask[snr_radio > snr_limits_radio[TELESCOPE]] = True
        print(f"  Mask fraction after ps : {np.sum(mask)/mask.size:.3f}")
        # write out
        fname_out_hdf5 = os.path.join(
            outdir, f"mask_{complexity}complexity.{telescope}.f{freq:03}.h5"
        )
        write_healpix(
            fname_out_hdf5, mask, dtype=np.int16, coord="C", nest=True, overwrite=True
        )
        print(f"Wrote {fname_out_hdf5}")
        fname_out_fits = os.path.join(
            outdir, f"mask_{complexity}complexity.{telescope}.f{freq:03}.fits"
        )
        hp.write_map(
            fname_out_fits, mask, dtype=np.int16, coord="C", nest=True, overwrite=True
        )
        print(f"Wrote {fname_out_fits}", flush=True)
