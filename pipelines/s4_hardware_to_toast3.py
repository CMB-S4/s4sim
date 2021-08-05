#!/usr/bin/env/python

# Copyright (c) 2020-2020 CMB-S4 Collaboration.
# Full license can be found in the top level "LICENSE" file.

"""Translate an S4 hardware map on disk into an HDF5 supported by TOAST3
"""

from collections import OrderedDict
import sys
import argparse

from astropy.table import QTable, Column
import astropy.units as u
import numpy as np
import toast.qarray as qa
from toast.instrument import Focalplane

import s4sim.hardware as hardware


XAXIS, YAXIS, ZAXIS = np.eye(3)


def main():
    parser = argparse.ArgumentParser(
        description="This program simulates or reads a hardware model "
        "and writes out a TOAST3 focalpane file.",
        usage="s4_hardware_to_toast3 [options] (use --help for details)",
    )

    parser.add_argument(
        "--hardware",
        required=False,
        help="Input hardware file.  Provide --telescope if not providing a file.",
    )
    
    parser.add_argument(
        "--telescope",
        required=True,
        help="Input telescope (LAT0,...,LAT2, or SAT0,...,SAT5).  "
        "Not used if --hardware is provided.",
    )

    parser.add_argument(
        "--out",
        required=False,
        default="focalplane",
        help="Name (without extensions) of the output hardware file",
    )

    parser.add_argument(
        "--fsample",
        required=True,
        default=None,
        type=float,
        help="Sampling rate for the focalplane",
    )

    args = parser.parse_args()

    if args.hardware is None:
        if args.telescope is None:
            raise RuntimeError("Must select a telescope if not providing a hardware file.")
        print(f"Simulating hardware for {args.telescope}...", end="", flush=True)
        hw = hardware.get_example()
        hw.data["detectors"] = hardware.sim_telescope_detectors(hw, args.telescope)
    else:
        if args.telescope is not None:
            raise RuntimeError("Must not select a telescope if providing a hardware file.")
        print(f"Loading hardware from {args.hardware}...", end="", flush=True)
        hw = hardware.Hardware(args.hardware)
    print(" Done!")

    # Each band must be written to a separate focaplane file
    bands = {}
    for det_name, det_data in hw.data["detectors"].items():
        band = det_data["band"]
        """
        # Optional code for getting polarization sentisive angle
        quat = qa.norm(det_data["quat"])
        wx, wy, wz = qa.rotate(quat, ZAXIS)
        wdir = np.array([wx, wy, wz])
        posrot = qa.norm(qa.from_vectors(ZAXIS, wdir))
        angrot = qa.norm(qa.mult(qa.inv(posrot), quat))
        psi = np.degrees(qa.to_angles(angrot)[2]) % 180
        det_data["psi_pol"] = psi
        """
        if band not in bands:
            bands[band] = {"name" : []}
        bands[band]["name"].append(det_name)
        # Keys are wafer, ID, pixel, fwhm, pol, card, channel, coax and bias
        for key, value in det_data.items():
            if key not in bands[band]:
                bands[band][key] = []
            bands[band][key].append(value)

    for band_name, det_data in bands.items():
        n_det = len(det_data["name"])

        band_data = hw.data["bands"][band_name]
        bandcenter = band_data["center"]
        bandwidth = band_data["high"] - band_data["low"]
        net = band_data["NET"] * 1e-6  # to K rts
        fknee = band_data["fknee"] * 1e-3  # to Hz
        fmin = band_data["fmin"] * 1e-3  # to Hz
        # The `alpha` in the hardware map includes atmosphere and
        # is much too large for instrumental noise
        # alpha = band_data["alpha"]
        alpha = 1.0
        A = band_data["A"]
        C = band_data["C"]
        pol_leakage = 0

        ones = np.ones(n_det)
        columns = [
            Column(name="name", data=det_data["name"]),
            Column(name="pol_leakage",data= ones * pol_leakage, unit=None),
            #Column(name="fwhm", data=ones * fwhm, unit=u.arcmin),
            Column(name="psd_fmin", data=ones * fmin, unit=u.Hz),
            Column(name="psd_fknee", data=ones * fknee, unit=u.Hz),
            Column(name="psd_alpha", data=ones * alpha, unit=None),
            Column(name="psd_net", data=ones * net, unit=(u.K * u.second ** .5)),
            Column(name="bandcenter", data=ones * bandcenter, unit=u.GHz),
            Column(name="bandwidth", data=ones * bandwidth, unit=u.GHz),
            Column(name="A", data=ones * A, unit=None),
            Column(name="C", data=ones * C, unit=None),
        ]
        for key, value in det_data.items():
            unit = None
            if key == "name":
                continue
            if key == "fwhm":
                unit = u.arcmin
            columns.append(Column(name=key, data=value, unit=unit))

        det_table = QTable(columns)

        fp = Focalplane(detector_data=det_table, sample_rate=args.fsample * u.Hz)
        fname = f"{args.out}_{args.telescope}_{band_name}.h5"
        fp.write(fname)
        print(f"Wrote {fname}")
    
    return


if __name__ == "__main__":
    main()
