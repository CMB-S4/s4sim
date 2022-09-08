from glob import glob
from datetime import datetime
import os
import sys

import astropy.units as u
import h5py
import matplotlib.pyplot as plt
import numpy as np

from toast.weather import SimWeather

telescope = "LAT0_CHLAT"
band = "f030"

rootdir = f"/global/cfs/cdirs/cmbs4/dc/dc1/staging/noise_sim/outputs"


all_pwv = []
all_time = []
all_rms_raw = []
all_rms_filtered = []

nrow, ncol = 1, 2
iplot = 0
fig = plt.figure(figsize=[6 * ncol, 6 * nrow])
# iplot += 1
# ax1 = fig.add_subplot(nrow, ncol, iplot)
# iplot += 1
# ax2 = fig.add_subplot(nrow, ncol, iplot)
weather = None
failed = []
for indir in glob(f"{rootdir}/{telescope}/{band}/{telescope}_split_schedule_????"):
    try:
        obsname = os.path.basename(glob(f"{indir}/mapmaker_*_hits.h5")[0])
    except:
        print(f"No hit map: {indir}")
        failed.append(indir)
        continue
    obsname = obsname.replace("mapmaker_", "").replace("_hits.h5", "")

    raw = h5py.File(f"{indir}/raw_statistics_{obsname}.h5", "r")
    filtered = h5py.File(f"{indir}/filtered_statistics_{obsname}.h5", "r")

    datadump = h5py.File(glob(f"{indir}/{obsname}_*.h5")[0], "r")
    attrib = datadump["instrument"].attrs
    max_pwv = attrib["site_weather_max_pwv"]
    if max_pwv == "NONE":
        max_pwv = None
    else:
        max_pwv *= u.mm
    site_uid = attrib["site_uid"]
    site_name = attrib["site_weather_name"]
    realization = attrib["site_weather_realization"]
    time = attrib["site_weather_time"]
    date = datetime.utcfromtimestamp(time)
    if weather is None:
        weather = SimWeather(
            time=date, site_uid=site_uid, realization=realization, name=site_name, max_pwv=max_pwv
        )
    else:
        weather.set(time=date, site_uid=site_uid, realization=realization)
    pwv = weather.pwv # .to_value(u.mm)
    

    good = raw["variance"][:] != 0
    idet = np.arange(good.size)

    """
    label = f"{pwv.to_value(u.mm):.1f} mm"
    ax1.plot(idet[good], raw["variance"][good], '.', label=label)
    ax2.plot(idet[good], filtered["variance"][good], '.', label=label)
    """

    all_time.append(time)
    all_pwv.append(pwv.to_value(u.mm))
    all_rms_raw.append(np.median(np.sqrt(raw["variance"][good])))
    all_rms_filtered.append(np.median(np.sqrt(filtered["variance"][good])))

#ax1.legend(loc="best")

iplot += 1
ax = fig.add_subplot(nrow, ncol, iplot)
ax.plot(all_pwv, all_rms_raw, '.')
ax.plot(all_pwv, all_rms_filtered, '.')
ax.set_xlabel("PWV [mm]")
ax.set_ylabel("RMS [K$_\mathrm{CMB}$]")
ax.set_title("Raw RMS")

iplot += 1
ax = fig.add_subplot(nrow, ncol, iplot)
ax.plot(all_pwv, all_rms_filtered, '.')
ax.set_xlabel("PWV [mm]")
ax.set_ylabel("RMS [K$_\mathrm{CMB}$]")
ax.set_title("Filtered RMS")

fig.savefig(f"rms_vs_pwv_{telescope}_{band}.png")

plt.show()
