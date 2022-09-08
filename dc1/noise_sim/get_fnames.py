#!/usr/bin/env python

# Return the list of observing schedules that do not have existing log files attached

import glob
import os
import sys

schedule_dir = sys.argv[1]
log_dir = sys.argv[2]
TELESCOPE = sys.argv[3]
bands = sys.argv[4:]

pattern = f"{schedule_dir}/split_schedule_????.txt"
all_schedules = sorted(glob.glob(pattern))
nschedule = len(all_schedules)

starts = []
for band in bands:
    rootname = TELESCOPE + "_" + os.path.basename(all_schedules[0]).replace(".txt", "")
    logfile = f"{log_dir}/{band}/{rootname}_{band}.log"
    if not os.path.isfile(logfile):
        # No log files exist
        starts.append(0)
        continue

    rootname = TELESCOPE + "_" + os.path.basename(all_schedules[nschedule - 1]).replace(".txt", "")
    logfile = f"{log_dir}/{band}/{rootname}_{band}.log"
    if os.path.isfile(logfile):
        # all log files exist
        starts.append(nschedule)
        continue

    # Binary search for the limit

    lower = 0
    upper = nschedule - 1

    while True:
        middle = int(lower + (upper - lower) / 2)
        rootname = TELESCOPE + "_" + os.path.basename(all_schedules[middle]).replace(".txt", "")
        logfile = f"{log_dir}/{band}/{rootname}_{band}.log"
        if os.path.isfile(logfile):
            lower = middle
        else:
            upper = middle
        if upper == lower + 1:
            starts.append(upper)
            break

schedules = all_schedules[min(starts):]
nschedule = len(schedules)

print(" ".join(schedules))
