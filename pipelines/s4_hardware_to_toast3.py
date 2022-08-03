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
import toast.io

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
        help="Input telescope (LAT0,...,LAT2, or SAT0,...,SAT5, or CHSAT0,...,CHSAT5 ).  "
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
        required=False,
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
        if args.fsample is not None:
            fsample = args.fsample
        else:
            # Hard-coded sampling rates
            if args.telescope.startswith("SAT") or args.telescope.startswith("CHSAT"):
                fsample = 20
            elif band_name.endswith("f220") or band_name.endswith("f280"):
                fsample = 440
            else:
                fsample = 220
                
        n_det = len(det_data["name"])
        print(
            f"\nFound {n_det} detectors in band {band_name}, sampling at {fsample} Hz"
        )

        tubes = []
        wafer_to_tube = {}
        wafer_to_FOV_cut = {}
        for wafer in set(det_data["wafer"]):
            # Determine which tube has this wafer
            for tube_name, tube_data in hw.data["tubes"].items():
                if wafer in tube_data["wafers"]:
                    break
            else:
                raise RuntimeError(f"Could not match wafer {wafer} to a tube.")
            wafer_to_tube[wafer] = tube_name
            tube_data = hw.data["tubes"][tube_name]
            if "FOV_cut" in tube_data:
                wafer_to_FOV_cut[wafer] = np.radians(tube_data["FOV_cut"])

        # Create the tube column
        tubes = [wafer_to_tube[wafer] for wafer in det_data["wafer"]]
        det_data["tube"] = tubes

        # Make pixel ID unique
        for idet, (tube, wafer, pixel) in enumerate(
                zip(det_data["tube"], det_data["wafer"], det_data["pixel"])
        ):
            det_data["pixel"][idet] = f"{tube}-{wafer}-{pixel}"

        # Cut detectors based on FOV_cut
        cut = np.ones(n_det, dtype=bool)
        for idet, (wafer, quat) in enumerate(
                zip(det_data["wafer"], det_data["quat"])
        ):
            if wafer in wafer_to_FOV_cut:
                theta_max = wafer_to_FOV_cut[wafer] / 2
                theta, phi = qa.to_position(quat)
                if theta > theta_max:
                    cut[idet] = False
        print(f"INFO: FOV_cut rejects {(1 - np.sum(cut) / n_det) * 100:.2f} % pixels")
        print(
            f"INFO: FOV_cut: {cut.size:8} - {cut.size - np.sum(cut):8}"
            f" = {np.sum(cut):8} detectors"
        )
        n_det = np.sum(cut)

        band_data = hw.data["bands"][band_name]
        bandcenter = band_data["center"]
        bandwidth = band_data["high"] - band_data["low"]
        net = band_data["NET"] * 1e-6  # to K rts
        if "NET_corr" in band_data:
            net_corr = band_data["NET_corr"]
            print(f"Scaling single detector NETs with the correlation factor: {net_corr}")
            net *= net_corr
        fknee = band_data["fknee"] * 1e-3  # to Hz
        fmin = band_data["fmin"] * 1e-3  # to Hz
        # The `alpha` in the hardware map includes atmosphere and
        # is much too large for instrumental noise
        # alpha = band_data["alpha"]
        alpha = 1.0
        A = band_data["A"]
        C = band_data["C"]
        pwv_poly = band_data["pwv_poly"]
        pol_leakage = 0

        ones = np.ones(n_det)
        names = []
        for name, flag in zip(det_data["name"], cut):
            if flag:
                names.append(name)
        columns = [
            Column(name="name", data=names),
            Column(name="pol_leakage", data= ones * pol_leakage, unit=None),
            #Column(name="fwhm", data=ones * fwhm, unit=u.arcmin),
            Column(name="psd_fmin", data=ones * fmin, unit=u.Hz),
            Column(name="psd_fknee", data=ones * fknee, unit=u.Hz),
            Column(name="psd_alpha", data=ones * alpha, unit=None),
            Column(name="psd_net", data=ones * net, unit=(u.K * u.second ** .5)),
            Column(name="bandcenter", data=ones * bandcenter, unit=u.GHz),
            Column(name="bandwidth", data=ones * bandwidth, unit=u.GHz),
            Column(name="elevation_noise_a", data=ones * A, unit=None),
            Column(name="elevation_noise_c", data=ones * C, unit=None),
            Column(name="pwv_noise_a0", data=ones * pwv_poly[0], unit=None),
            Column(name="pwv_noise_a1", data=ones * pwv_poly[1], unit=None),
            Column(name="pwv_noise_a2", data=ones * pwv_poly[2], unit=None),
        ]

        for key, value in det_data.items():
            unit = None
            if key == "name":
                continue
            if key == "fwhm":
                unit = u.arcmin
            columns.append(Column(name=key, data=np.array(value)[cut], unit=unit))

        det_table = QTable(columns)

        fp = Focalplane(detector_data=det_table, sample_rate=fsample * u.Hz)
        fname = f"{args.out}_{args.telescope}_{band_name}.h5"
        with toast.io.H5File(fname, "w", comm=None, force_serial=True) as f:
            fp.save_hdf5(f.handle)
        print(f"Wrote {fname}")

    return


if __name__ == "__main__":
    main()
