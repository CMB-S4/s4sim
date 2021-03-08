import os
import sys

import healpy as hp
import numpy as np

import s4sim.hardware as hardware


print("Simulating hardware", flush=True)
hw = hardware.get_example()

print("Simulating detectors", flush=True)
hw.data["detectors"] = hardware.sim_telescope_detectors(hw, "SAT0")

# sigma = net * sqrt(fsample) => net = sigma / sqrt(fsample)
net_in = np.sqrt(1 / 1e-3) / np.sqrt(20)

for pixel_type in "LF", "MFL", "MFH", "HF":
    for iband in 1, 2:
        band = f"{pixel_type}S{iband}"
        fname_in = os.path.join(
            "out/00000000/temp",
            f"pole_noise_SAT_{band}_telescope_all_time_all_wcov.fits.gz"
        )
        fname_out = os.path.join(
            "out/00000000/",
            f"pole_noise_SAT_{band}_telescope_all_time_all_wcov.fits.gz"
        )
        bdata = hw.data["bands"][band]
        net = bdata["NET"] * 1e-6  # in K_CMB
        sigma = net * np.sqrt(20)
        print(fname_in, net)
        print("Loading", fname_in)
        w = hp.read_map(fname_in, None, nest=True)
        w *= sigma ** 2 / 1e-3
        print("Writing", fname_out)
        hp.write_map(fname_out, w, nest=True, overwrite=True)
