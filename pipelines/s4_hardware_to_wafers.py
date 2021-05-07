#!/usr/bin/env/python

# Copyright (c) 2020-2020 CMB-S4 Collaboration.
# Full license can be found in the top level "LICENSE" file.
"""Trim a hardware model to include only some detectors.
"""

from collections import OrderedDict
import sys
import argparse

import numpy as np
import toast.qarray as qa

from s4sim.hardware import Hardware


XAXIS, YAXIS, ZAXIS = np.eye(3)


def main():
    parser = argparse.ArgumentParser(
        description="This program reads a hardware model from disk "
        "and writes out per-wafer geometry.",
        usage="s4_hardware_trim [options] (use --help for details)",
    )

    parser.add_argument(
        "--hardware", required=True, default=None, help="Input hardware file"
    )

    parser.add_argument(
        "--out",
        required=False,
        default="trimmed",
        help="Name (without extensions) of the output hardware file",
    )

    args = parser.parse_args()

    print("Loading hardware from {}...".format(args.hardware), flush=True)
    hw = Hardware(args.hardware)

    wafers = {}

    for det_name, det_data in hw.data["detectors"].items():
        wafer = det_data["wafer"]
        pol = det_data["pol"]
        quat = qa.norm(det_data["quat"])
        if wafer not in wafers:
            wafers[wafer] = {}
        wx, wy, wz = qa.rotate(quat, ZAXIS)
        wdir = np.array([wx, wy, wz])
        posrot = qa.norm(qa.from_vectors(ZAXIS, wdir))
        angrot = qa.norm(qa.mult(qa.inv(posrot), quat))
        psi = np.degrees(qa.to_angles(angrot)[2]) % 180
        wafers[wafer][det_name] = np.array([wx, wy, psi])

    for wafer_name in sorted(wafers):
        wafer_data = wafers[wafer_name]
        ndet = len(wafer_data)
        det_data = np.empty([ndet, 3])
        for idet, det_name in enumerate(wafer_data):
            det_data[idet] = wafer_data[det_name]
        xoffset, yoffset, psioffset = np.mean(det_data, 0)
        fname = f"wafer_{wafer_name}.txt"
        with open(fname, "w") as fout:
            for det_name in sorted(wafer_data):
                x, y, psi = wafer_data[det_name]
                fout.write(f"{det_name} {x - xoffset} {y - yoffset} {psi}\n")

    return


if __name__ == "__main__":
    main()
