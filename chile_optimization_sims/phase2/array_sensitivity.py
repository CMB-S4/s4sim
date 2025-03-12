import glob
import os
import sys

import astropy.units as u
import numpy as np
import toast


yield_ = 0.8

print("All NETs given for polarization, not temperature")


for band in 30, 40, 85, 90, 95, 145, 150, 155, 220, 280:
    fname = glob.glob(f"focalplanes/focalplane_SAT?_SAT_f{band:03}_ST?.h5")[0]
    if band in [30, 40]:
        ntube = 3
    elif band in [85, 95, 145, 155]:
        ntube = 9
    elif band in [90, 150]:
        ntube = 18
    elif band in [220, 280]:
        ntube = 6
    else:
        raise RuntimeError(f"Unknown band: {band}")
    focalplane = toast.instrument.Focalplane()
    with toast.io.H5File(fname, "r") as f:
        focalplane.load_hdf5(f.handle)
    net = focalplane.detector_data["psd_net"].to_value(u.uK*u.s**.5)
    net *= np.sqrt(2) # T->P
    ndet = net.size * ntube
    net = np.mean(net)
    net_array = net / np.sqrt(ndet * yield_)
    print(f"SAT {band:03} : NET(array) = {net_array:8.3f} uKrts NET(det) = {net:8.3f} ndet = {ndet}")


for band in 20, 30, 40, 90, 150, 220, 280:
    fname = f"focalplanes/focalplane_LAT0_CHLAT_f{band:03}.h5"
    focalplane = toast.instrument.Focalplane()
    with toast.io.H5File(fname, "r") as f:
        focalplane.load_hdf5(f.handle)
    net = focalplane.detector_data["psd_net"].to_value(u.uK*u.s**.5)
    net *= np.sqrt(2) # T->P
    ndet = net.size
    net = np.mean(net)
    net_array = net / np.sqrt(ndet * yield_)
    print(f"LAT {band:03} : NET(array) = {net_array:8.3f} uKrts NET(det) = {net:8.3f} ndet = {ndet}")
