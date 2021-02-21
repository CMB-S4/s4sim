#!/usr/bin/env python

# Copyright (c) 2020-2020 CMB-S4 Collaboration..
# Full license can be found in the top level "LICENSE" file.

# This script will convert TOML hardware maps into Python pickle
# format for efficient access

import os
import pickle
import sys

from s4sim import hardware


if len(sys.argv) < 2:
    print("Usage: pickle_hardware.py <TOML file1> [<TOML file2>] ...")
else:
    for fname_in in sys.argv[1:]:
        print(f"Loading {fname_in}")
        hw = hardware.Hardware(fname_in)
        
        fname_out = fname_in.replace(".toml", "").replace(".gz", "") + ".pkl"
        print(f"Writing {fname_out}")
        with open(fname_out, "wb") as fout:
            pickle.dump(hw, fout)
