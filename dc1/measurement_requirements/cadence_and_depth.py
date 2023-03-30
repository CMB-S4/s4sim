import os
import pickle
import sys

import astropy.units as u
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import c, h, k

from toast.pixels_io_healpix import read_healpix


fname_out = "cadence_and_depth.pck"
if os.path.isfile(fname_out):
    with open(fname_out, "rb") as f:
        common, depths = pickle.load(f)
    """
    # TEMPORARY
    temp = common.copy()
    common = {}
    for date, fsky in zip(temp["date"], temp["fsky"]):
        common[date] = fsky
    temp = depths.copy()
    depths = {}
    for date, depth in zip(temp["date"], temp["depth"]):
        depths[date] = depth
    """
else:
    common = {}
    depths = {}

nrow, ncol = 1, 2
fig = plt.figure(figsize=[4 * ncol, 3 * nrow])

x, y = [], []
for key, value in depths.items():
    x.append(key)
    y.append(value)
y = np.array(y)
ax = fig.add_subplot(nrow, ncol, 1)
doy = np.arange(y.size) + 1
ax.plot(doy, y, '.', label="Best 0.25 fsky")
ax.set_xlabel("DOY")
ax.set_ylabel("Depth [mJy]")
ax.axhline(10.0, color="k", linestyle="--", label="MR 4.1")
ax.set_title("Daily depth")
ax.legend(loc="best")
ax.set_ylim([0, 20])

x, y = [], []
for key, value in common.items():
    x.append(key)
    y.append(value)
y = np.array(y)
ax = fig.add_subplot(nrow, ncol, 2)
doy = np.arange(y.size) + 1
ax.plot(doy, y, '.', label="Last 5 days")
ax.set_xlabel("DOY")
ax.set_ylabel("fsky")
ax.axhline(0.25, color="k", linestyle="--", label="MR 4.1")
ax.set_title("5-day common sky fraction")
ax.legend(loc="best")

fig.tight_layout()
plt.savefig("cadence_and_depth.pdf")

TCMB = 2.72548
freq = 150 * u.GHz
fwhm = 1.4 * u.arcmin
sigma = fwhm / np.sqrt(8 * np.log(2))
solid_angle = 2 * np.pi * sigma**2

#fname_schedule = "../scan_strategy/chile_lat/schedules/chile_schedule_lat.pruned.txt"
fname_schedule = "../scan_strategy/chile_lat/schedules/chile_schedule_lat.txt"

# /global/cfs/cdirs/cmbs4/dc/dc1/staging/noise_sim/outputs_rk/LAT0_CHLAT/f150/RISING_SCAN_40-231-12/mapmaker_RISING_SCAN_40-231-12_hits.h5
#indir = "/global/cfs/cdirs/cmbs4/dc/dc1/staging/noise_sim/outputs_rk/LAT0_CHLAT/f150"
indir = "/global/cfs/cdirs/cmbs4/dc/dc1/staging/foreground_sim/outputs_rk/LAT0_CHLAT/f150"


schedule = []
dates = set()
with open(fname_schedule, "r") as file_in:
    for iline, line in enumerate(file_in):
        if iline < 3:
            continue
        parts = line.split()
        start_date = parts[0]
        dates.add(start_date)
        start_time = parts[1]
        start = f"{start_date} {start_time}"
        target = parts[7]
        scan = parts[21]
        subscan = parts[22]
        observation = f"{target}-{scan}-{subscan}"
        schedule.append((start, observation))

nside = 4096
npix = 12 * nside ** 2
limit = int(npix * 0.25)  # Measured over best 25% of the sky
pix_area = hp.nside2pixarea(nside)
year = 2027
n = 5
daily_invcovs = []
total_hits = np.zeros(npix, dtype=np.int32)

# Delta-function bandpass

nu_cmb = k * TCMB / h
alpha = 2 * k**3 * TCMB**2 / h**2 / c**2
x = freq.to_value(u.Hz) / nu_cmb
db_dt = alpha * x**4 * np.exp(x) / (np.exp(x) - 1) ** 2
kcmb2Jysr = db_dt * 1e26
kcmb2mJysr = kcmb2Jysr * 1e3

sys.exit()

for month in range(12, 13):
    for day in range(1, 32):
        date = f"{year}-{month:02}-{day:02}"
        if date not in dates:
            # not a scheduled date
            continue
        print(f"\n{date}", flush=True)
        if len(daily_invcovs) == n:
            oldest = daily_invcovs.pop(0)
            total_hits -= (oldest != 0)
        daily_invcov = np.zeros(npix)
        for start, observation in schedule:
            if start.split()[0] != date:
                continue
            fname_invcov = os.path.join(
                indir,
                observation,
                f"mapmaker_{observation}_invcov.h5"
            )
            if os.path.isfile(fname_invcov):
                print(f"Reading {fname_invcov}", flush=True)
                invcov = read_healpix(fname_invcov, [0], dtype=np.float32, nest=True)[0]
                daily_invcov += invcov
        good = daily_invcov != 0
        ngood = np.sum(good)
        if ngood < limit:
            depth = np.inf
        else:
            depth = np.sqrt(1 / daily_invcov[good]) * kcmb2mJysr \
                * solid_angle.to_value(u.steradian)
            depth = np.sort(depth)[limit]
            print(f"Depth = {depth:.2f} mJy", flush=True)
        depths["date"] = depth
        daily_invcovs.append(daily_invcov)
        total_hits += (daily_invcov != 0)
        if len(daily_invcovs) == n:
            fsky = np.sum(total_hits == n) / npix
            common[date] = fsky
            print(f"Common fsky {fsky:.3f}", flush=True)
        with open(fname_out, "wb") as f:
            pickle.dump([common, depths], f)
        print(f"Wrote {fname_out}", flush=True)
