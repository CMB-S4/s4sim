#!/usr/bin/env python

from glob import glob
import os
import sys

import numpy as np

for tele in "chlat", "splat", "spsat":
    print(f"\n{tele}\n")
    fnames = glob(f"split_schedules/{tele}/split_schedule_*.txt")
    fnames_done = glob(f"split_schedules/{tele}/done/split_schedule_*.txt")
    nsplit = len(fnames) + len(fnames_done)
    print(f"Found {nsplit} split schedule files")

    for freq in 30, 40, 90, 150, 220, 280:
        fnames = glob(f"logs/{tele}/f{freq:03}/*.log")
        nlog = len(fnames)
        print(f"{freq:3} GHz is {nlog:5} / {nsplit:5} = {nlog / nsplit * 100:6.3f}% complete")
