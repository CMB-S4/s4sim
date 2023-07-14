import os
import sys

import numpy as np
import healpy as hp
from toast.pixels_io_healpix import write_healpix


outdir = "input_maps"
os.makedirs(outdir, exist_ok=True)

bands = {
    "CHLAT" : ("chlat", 4096, [30, 40, 90, 150, 220, 280]),
    "SPLAT" : ("splat", 4096, [20, 30, 40, 90, 150, 220, 280]),
    "SAT" : ("spsat", 512, [30, 40, 85, 95, 145, 155, 220, 280]),
}

for TELESCOPE, (telescope, nside, freqs) in bands.items():
    for freq in freqs:
        fname_out = os.path.join(outdir, f"cmb_lensing.{telescope}.f{freq:03}.h5")
        if os.path.isfile(fname_out):
            print(f"Output file exists: {fname_out}")
            continue
        fname_in = f"/pscratch/sd/z/zonca/cmbs4/202305_dc0/combined_cmb_lensing_signal/cmbs4_combined_cmb_lensing_signal_uKCMB_{TELESCOPE}_f{freq:03}_nside{nside}.fits"
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
        write_healpix(fname_out, m * 1e-6, coord="C", nest=True)
        print(f"Wrote {fname_out}")
