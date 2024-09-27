import argparse

from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
import astropy.units as u
import numpy as np


parser = argparse.ArgumentParser()
parser.add_argument("ifile")
parser.add_argument("ofile")
args = parser.parse_args()

# config.set("iers_fallback", "nearest")
nsamp = 100

# Read in the file, skipping comments
with open(args.ifile, "r") as ifile:
    lines = [line for line in ifile if len(line) != 0]

header = lines[:3]
lines = lines[3:]

parts = header[1].split()
lat = float(parts[2]) * u.degree
lon = float(parts[3]) * u.degree
height = float(parts[4]) * u.m
loc = EarthLocation(lat=lat, lon=lon, height=height)

scan_ras = []
for line in lines:
    print(line.rstrip())
    toks = line.split()
    time1 = Time(f"{toks[0]} {toks[1]}", format="iso")
    time2 = Time(f"{toks[2]} {toks[3]}", format="iso")
    az1 = float(toks[6])  # degrees
    az2 = float(toks[7])  # degrees
    el = float(toks[8])  # degrees
    rising = az1 < 180
    azs = np.linspace(az1, az2, nsamp) * u.degree
    els = np.full(nsamp, el) * u.degree
    radec1 = SkyCoord(alt=els, az=azs, obstime=time1, frame="altaz", location=loc).icrs
    radec2 = SkyCoord(alt=els, az=azs, obstime=time2, frame="altaz", location=loc).icrs
    ras = np.concatenate([radec1.ra.to_value(u.deg), radec2.ra.to_value(u.deg)])
    try:
        ras = np.unwrap(ras, period=360)
    except TypeError:
        ras = np.unwrap(ras, discont=360)
    middle = np.mean([np.amin(ras), np.amax(ras)])
    if rising:
        scan_ras.append(middle)
    else:
        # Sort the setting scans separate from the rising ones
        scan_ras.append(-middle)

# Sort after accounting for any wrapping issues
order = np.argsort(scan_ras)

# Finally output the lines we read in in the new order
with open(args.ofile, "w") as ofile:
    for line in header:
        ofile.write(line)
    for i in order:
        # ofile.write(lines[i].rstrip() + f"{scan_ras[i]:10.3f}\n")
        ofile.write(lines[i])
