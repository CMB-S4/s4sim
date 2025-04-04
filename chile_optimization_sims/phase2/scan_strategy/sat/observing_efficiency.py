import os
import sys

import ephem
import astropy.units as u
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np

import toast
from toast import qarray as qa
from toast.coordinates import to_DJD


# angles = [90, 80, 70, 60, 50, 40]
angles = [90]
flavors = ["max", "max+bk"]
print(f"flavor angle all/year all/season south/year south/season")
for angle in angles:
    for flavor in flavors:
        tsouth = 0
        tsouth_full = 0
        tnorth = 0
        tnorth_full = 0
        schedule = toast.schedule.GroundSchedule()
        schedule.read(f"schedule_sat.sun{angle}{flavor}.txt")
        t0 = schedule.scans[0].start.timestamp()
        for scan in schedule.scans:
            tstart = scan.start.timestamp()
            tstop = scan.stop.timestamp()
            tdelta = tstop - tstart
            in_season = (tstart - t0) / 86400 > 90
            if "South" in scan.name:
                tsouth_full += tdelta
                if in_season:
                    tsouth += tdelta
            else:
                tnorth_full += tdelta
                if in_season:
                    tnorth += tdelta
        year = 86400 * 365
        season = 86400 * (365 - 90)
        print(f"{flavor:8} {angle:3}"
              f" {(tsouth_full + tnorth_full) / year:.3f} {(tsouth + tnorth) / season:.3f}"
              f" {tsouth_full / year:.3f} {tsouth / season:.3f}"
              )
