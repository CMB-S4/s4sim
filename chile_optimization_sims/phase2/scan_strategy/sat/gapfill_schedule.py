# Fill gaps in one schedule with entries from the other
# Output is a supplemental schedule that can be merged with the master schedule

import datetime
import os
import sys

import healpy as hp
import numpy as np


if len(sys.argv) < 4:
    print(f"Usage: {sys.argv[0]} <master_schedule> <supplementing_schedule> <output_schedule>")
    sys.exit()
    
fname_master = sys.argv[1]
fname_supplement = sys.argv[2]
fname_out = sys.argv[3]

# load the input schedules
header = []
master = []
with open(fname_master, "r") as f:
    for iline, line in enumerate(f):
        if iline < 3:
            header.append(line)
            continue
        parts = line.split()
        sstart = f"{parts[0]} {parts[1]}"
        sstop = f"{parts[2]} {parts[3]}"
        start = datetime.datetime.fromisoformat(sstart).timestamp()
        stop = datetime.datetime.fromisoformat(sstop).timestamp()
        master.append([start, stop, line])

supplement = []
with open(fname_supplement, "r") as f:
    for iline, line in enumerate(f):
        if iline < 3:
            continue
        parts = line.split()
        sstart = f"{parts[0]} {parts[1]}"
        sstop = f"{parts[2]} {parts[3]}"
        start = datetime.datetime.fromisoformat(sstart).timestamp()
        stop = datetime.datetime.fromisoformat(sstop).timestamp()
        supplement.append([start, stop, line])

with open(fname_out, "w") as f:
    for line in header:
        f.write(line)
    nline = len(master)
    for iline in range(nline - 1):
        # f.write(master[iline][2])
        gap_start = master[iline][1]
        gap_stop = master[iline + 1][0]
        if gap_stop - gap_start > 3600:
            # see if we can fill the gap from the other schedule
            for start, stop, line in supplement:
                if stop < gap_start:
                    continue
                if start > gap_stop:
                    break
                if gap_start < start and stop < gap_stop:
                    f.write(line)
    # f.write(master[-1][2])
