#!/usr/bin/env python

from glob import glob
import os
import sys

import numpy as np

tele = "chlat"
fnames = glob(f"split_schedules/{tele}/split_schedule_*.txt")
nsplit = len(fnames)
print(f"Found {nsplit} split schedule files")

for freq in 30, 40, 90, 150, 220, 280:
    fnames = glob(f"logs/chlat/f{freq:03}/*.log")
    nlog = len(fnames)
    print(f"{freq:3} GHz is {nlog:5} / {nsplit:5} = {nlog / nsplit * 100:6.3f}% complete")
