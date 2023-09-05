# This script parses an observing schedule and identifies the
# geometrically different scans

import argparse
import os
import sys

import astropy.units as u
import numpy as np

import toast
from toast.coordinates import to_UTC
import toast.schedule


# Argument parsing

parser = argparse.ArgumentParser(description="Find the Irreducible Scanset")

parser.add_argument(
    "schedule", nargs="+", help="One or more schedules"
)

parser.add_argument(
    "--pole_mode",
    required=False,
    default=False,
    action="store_true",
    help="In Pole mode rising and setting scans are equivalent",
)

parser.add_argument(
    "--asymmetric_focalplane",
    required=False,
    default=False,
    action="store_true",
    help="Treat boresight angles separated by 180 degrees as different",
)

parser.add_argument(
    "--tol",
    default=0.1,
    type=float,
    help="Tolerance in elevation and boresight angle",
)

args = parser.parse_args()

# Initialize

comm, procs, rank = toast.get_world()

# Load the schedule file

schedule = toast.schedule.GroundSchedule()
schedule.read(args.schedule[0], comm=comm)

irreducible = {}
class ScanData():
    def __init__(self):
        self.count = 0
        self.times = []
        self.azranges = []
        self.scans = []
        self.subscans = []

    def __repr__(self):
        return f"Count = {self.count}, times = {self.times[0][0]} -- {self.times[-1][1]}"

ntot = 0
subscans = None
for scan in schedule.scans:
    if int(scan.subscan_indx) != 0 and not args.pole_mode:
        # Only count the first subscans and assume that whole scans are being scheduled
        subscans.append(scan)
        continue
    name = scan.name
    el = scan.el.to_value(u.deg)
    angle = scan.boresight_angle.to_value(u.deg)
    if not args.asymmetric_focalplane:
        angle %= 180
    # Apply the tolerance
    el = int(el / args.tol) * args.tol
    angle = int(angle / args.tol) * args.tol
    if args.pole_mode:
        rising = True
    else:
        rising = scan.rising
    if name not in irreducible:
        irreducible[name] = {True : {}, False : {}}
    if el not in irreducible[name][rising]:
        irreducible[name][rising][el] = {}
    if angle not in irreducible[name][rising][el]:
        irreducible[name][rising][el][angle] = ScanData()
    irreducible[name][rising][el][angle].count += 1
    irreducible[name][rising][el][angle].times.append((scan.start, scan.stop))
    irreducible[name][rising][el][angle].azranges.append((scan.az_min, scan.az_max))
    irreducible[name][rising][el][angle].scans.append(scan)
    subscans = []
    irreducible[name][rising][el][angle].subscans.append(subscans)
    ntot += 1

# Print the results

ntot = 0
for target in irreducible:
    print(f"{target}")
    ntot_target = 0
    for rs in irreducible[target]:
        elevations = irreducible[target][rs]
        if len(elevations) == 0:
            continue
        print(f"  rising/setting = {rs}")
        els = sorted(elevations.keys())
        for el in els:
            print(f"    {el:8.3f}")
            angles = sorted(elevations[el].keys())
            for angle in angles:
                ntot_target += 1
                print(f"        {angle:8.3f} : {elevations[el][angle].count:8}")
    print(f"Irreducible scans {target} = {ntot_target}")
    ntot += ntot_target
print(f"Irreducible scans = {ntot}")


# Write a new schedule with the irreducible scans

fname_out = "irreducible_schedule.txt"
with open(fname_out, "w") as fout:
    site = schedule.site_name
    telescope = "TELESCOPE"
    lat = schedule.site_lat.to_value(u.deg)
    lon = schedule.site_lon.to_value(u.deg)
    alt = schedule.site_alt.to_value(u.m)
    fout.write(f"#Site Telescope Latitude[deg] Longitude[deg] Elevation[m]\n")
    fout.write(f"{site} {telescope} {lat} {lon} {alt}\n")

    # Concise schedule format
    fout_fmt0 = "#{:>20} {:>20} {:>8} {:35} {:>8} {:>8} {:>8} {:>5} {:>3}\n"
    fout_fmt = " {:>20} {:>20} {:8.2f} {:35} {:8.2f} {:8.2f} {:8.2f} {:5} {:3}\n"
    fout.write(
        fout_fmt0.format(
            "Start time UTC",
            "Stop time UTC",
            "Rotation",
            "Patch name",
            "Az min",
            "Az max",
            "El",
            "Pass",
            "Sub",
        )
    )

    for target in irreducible:
        for rs in irreducible[target]:
            elevations = irreducible[target][rs]
            if len(elevations) == 0:
                continue
            els = sorted(elevations.keys())
            for el in els:
                angles = sorted(elevations[el].keys())
                for angle in angles:
                    scandata = irreducible[target][rs][el][angle]
                    ndata = len(scandata.times)
                    mid = ndata // 2
                    start, stop = scandata.times[mid]
                    start = to_UTC(start.timestamp())
                    stop = to_UTC(stop.timestamp())
                    azmin, azmax = scandata.azranges[mid]
                    azmin = azmin.to_value(u.deg)
                    azmax = azmax.to_value(u.deg)
                    entry = fout_fmt.format(
                        start, stop, angle, target, azmin, azmax, el, ndata, 0
                    )
                    print(entry, end="")
                    fout.write(entry)
                    # Did this scan have subscans?
                    for subscan in scandata.subscans[mid]:
                        start = to_UTC(subscan.start.timestamp())
                        stop = to_UTC(subscan.stop.timestamp())
                        azmin = subscan.az_min.to_value(u.deg)
                        azmax = subscan.az_max.to_value(u.deg)
                        ind = int(subscan.subscan_indx)
                        entry = fout_fmt.format(
                            start, stop, angle, target, azmin, azmax, el, ndata, ind
                        )
                        print(entry, end="")
                        fout.write(entry)
print(f"Wrote {fname_out}")
