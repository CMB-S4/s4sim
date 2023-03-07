import os
import sys

import numpy as np
from healpy import read_map, smoothing
from toast.pixels_io_healpix import write_healpix


outdir = "input_maps"
os.makedirs(outdir, exist_ok=True)
mask_frac = 0.15
mask_fwhm = np.radians(1)
mask_lmax = 512

# CHLAT

fname_dipo = (
    "/global/cfs/cdirs/cmbs4/dm/mbs/"
    "202211_LAT_fg_cmb/4096/"
    "cmbs4_cmb_solar_dipole_uKCMB_nside4096.fits"
)
dipo = read_map(fname_dipo, nest=True)

for suffix in "", "_lowcomplexity", "_highcomplexity":
    for band, freq in [
            ("LFL1", 30),
            ("LFL2", 40),
            ("MFL1", 90),
            ("MFL2", 150),
            ("HFL1", 220),
            ("HFL2", 280),
    ]:
        fname_out = os.path.join(outdir, f"foreground{suffix}.chlat.f{freq:03}.h5")
        #if os.path.isfile(fname_out):
        #    print(f"Output file exists: {fname_out}")
        #    continue
        fname_in = (
            f"/global/cfs/cdirs/cmbs4/dm/mbs/"
            f"202211_LAT_fg_cmb/4096/combined_foregrounds_cmb_fixdip{suffix}/0000/"
            f"cmbs4_combined_foregrounds_cmb_fixdip{suffix}_uKCMB_LAT-{band}_nside4096_0000.fits"
        )
        if not os.path.isfile(fname_in):
            raise RuntimeError(f"Input file does not exist: {fname_in}")
        print(f"Reading {fname_in}")
        m = read_map(
            fname_in,
            None,
            nest=True,
            verbose=False,
            dtype=np.float32,
        )
        if not os.path.isfile(fname_out):
            write_healpix(fname_out, m * 1e-6, coord="C", nest=True)
            print(f"Wrote {fname_out}")

        # Processing mask
        # Subtract CMB
        fname_cmb = (
            f"/global/cfs/cdirs/cmbs4/dm/mbs/"
            f"202211_LAT_fg_cmb/4096/cmb/0000/"
            f"cmbs4_cmb_uKCMB_LAT-{band}_nside4096_0000.fits"
        )
        cmb = read_map(fname_cmb, nest=True)
        fg = m[0] - dipo - cmb
        fg_sorted = np.sort(fg)
        npix = fg.size
        lim = fg_sorted[int(npix * (1 - mask_frac))]
        mask = fg > lim
        mask = smoothing(mask, lmax=mask_lmax, fwhm=mask_fwhm, nest=True) > 0.25
        # mask radio sources separately
        fname_radio = (
            f"/global/cfs/cdirs/cmbs4/dm/mbs/"
            f"202211_LAT_fg_cmb/4096/radio/0000/"
            f"cmbs4_radio_uKCMB_LAT-{band}_nside4096_0000.fits"
        )
        radio = read_map(fname_radio, nest=True)
        radio_sorted = np.sort(radio)
        lim = radio_sorted[int(npix * .99)]
        mask[radio > lim] = True
        # write out
        fname_out = os.path.join(outdir, f"mask{int(mask_frac * 100):02}{suffix}.chlat.f{freq:03}.h5")
        write_healpix(fname_out, mask, dtype=np.int16, coord="C", nest=True, overwrite=True)
        print(f"Wrote {fname_out}")

    sys.exit()

# SPLAT

for band, freq in [
        ("ULFPL1", 20),
        ("LFPL1", 30),
        ("LFPL2", 40),
        ("MFPL1", 90),
        ("MFPL2", 150),
        ("HFPL1", 220),
        ("HFPL2", 280),
]:
    fname_out = os.path.join(outdir, f"cmb.splat.f{freq:03}.h5")
    if os.path.isfile(fname_out):
        print(f"Output file exists: {fname_out}")
        continue
    fname_in = (
        f"/global/cfs/cdirs/cmbs4/dm/dstool_202102/input_pysm"
        f"/4096/cmb_unlensed_solardipole/0000/"
        f"cmbs4_cmb_unlensed_solardipole_uKCMB_LAT-{band}_nside4096_0000.fits"
    )
    if not os.path.isfile(fname_in):
        raise RuntimeError(f"Input file does not exist: {fname_in}")
    print(f"Reading {fname_in}")
    m = read_map(
        fname_in,
        None,
        nest=True,
        verbose=False,
        dtype=np.float32,
    )
    write_healpix(fname_out, m * 1e-6, coord="C", nest=True)
    print(f"Wrote {fname_out}")

# SPSAT

for band, freq in [
        ("LFS1", 30),
        ("LFS2", 40),
        ("MFLS1", 85),
        ("MFLS2", 145),
        ("MFHS1", 95),
        ("MFHS2", 155),
        ("HFS1", 220),
        ("HFS2", 280),
]:
    fname_out = os.path.join(outdir, f"cmb.spsat.f{freq:03}.h5")
    if os.path.isfile(fname_out):
        print(f"Output file exists: {fname_out}")
        continue
    fname_in = (
        f"/global/cfs/cdirs/cmbs4/dm/dstool_202102/input_pysm"
        f"/512/cmb_unlensed_solardipole/0000/"
        f"cmbs4_cmb_unlensed_solardipole_uKCMB_SAT-{band}_nside512_0000.fits"
    )
    if not os.path.isfile(fname_in):
        raise RuntimeError(f"Input file does not exist: {fname_in}")
    print(f"Reading {fname_in}")
    m = read_map(
        fname_in,
        None,
        nest=True,
        verbose=False,
        dtype=np.float32,
    )
    write_healpix(fname_out, m * 1e-6, coord="C", nest=True)
    print(f"Wrote {fname_out}")

"""
python -c '
import numpy as np
from healpy import read_map;
from toast.pixels_io import write_healpix
m = read_map(
    None,
    nest=True,
    verbose=False,
    dtype=np.float32,
);
write_healpix(
    "cmb.spsat.f085.h5",
    m * 1e-6,
    coord="G",
    nest=True,
)
'
"""
