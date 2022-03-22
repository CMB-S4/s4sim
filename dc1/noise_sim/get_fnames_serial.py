#!/usr/bin/env python

# Return the list of observing schedules that do not have existing log files attached

import glob
import os
import sys

schedule_dir = sys.argv[1]
log_dir = sys.argv[2]
TELESCOPE = sys.argv[3]
bands = sys.argv[4:]

pattern = f"{schedule_dir}/*.txt"
all_schedules = sorted(glob.glob(pattern))
nschedule = len(all_schedules)

schedules = set()
for band in bands:
    for schedule in all_schedules:
        obs = os.path.basename(schedule).replace(".txt", "")
        logfile = f"{log_dir}/{band}/{obs}.log"
        if not os.path.isfile(logfile):
            schedules.add(schedule)

schedules = list(schedules)

print(" ".join(schedules))
