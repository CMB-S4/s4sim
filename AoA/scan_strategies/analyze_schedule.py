from glob import glob
import os
import sys

import matplotlib.pyplot as plt
import numpy as np


fig = plt.figure(figsize=[18, 12])
nrow, ncol = 2, 2
ax1 = fig.add_subplot(nrow, ncol, 1)
ax2 = fig.add_subplot(nrow, ncol, 2)
ax3 = fig.add_subplot(nrow, ncol, 3)

colors = ["tab:blue", "tab:orange", "tab:green", "tab:pink"]

if len(sys.argv) == 1:
    print("Usage: {sys.argv[0]} <input directory>")
indir = sys.argv[1]
fnames = sorted(glob(indir + "/*txt"))

for fname, color in zip(fnames, colors):
    print("\n{}".format(fname))
    schedule = np.genfromtxt(fname, skip_header=3).T

    arr = np.genfromtxt(
        fname,
        skip_header=3,
        usecols=[4, 5, 7, 8, 9, 10],
        dtype="f8,f8,S30,f8,f8,f8",
        names=["starts", "stops", "names", "az_mins", "az_maxs", "elevations"],
    )
    lengths = arr["stops"] - arr["starts"]
    throws = arr["az_maxs"] - arr["az_mins"]

    last_el = None
    for start, stop, az_min, az_max, el, name in zip(
        arr["starts"], arr["stops"], arr["az_mins"], arr["az_maxs"], arr["elevations"], arr["names"],
    ):
        if last_el is None:
            label = os.path.basename(fname).replace(".txt", "")
        else:
            label = None
        ax1.plot([start, stop], [el, el], color=color, label=label)
        if last_el is not None:
            ax2.plot([last_stop, start], [last_el, el], color=color, label=label)
        ax3.plot([az_min, az_max], [el, el], color=color, label=label)
        last_stop = stop
        last_el = el

    is_horizontal = np.array([
        name.decode().startswith("RISING_SCAN") or name.decode().startswith("SETTING_SCAN") for name in arr["names"]
    ])

    average_el = np.sum(arr["elevations"] * lengths) / np.sum(lengths)

    integration_time = np.sum(arr["stops"] - arr["starts"])
    integration_time_horizontal = np.sum(arr["stops"][is_horizontal] - arr["starts"][is_horizontal])
    schedule_time = arr["stops"][-1] - arr["starts"][0]
    t_steps = arr["starts"][1:] - arr["stops"][:-1]
    el_steps = np.diff(arr["elevations"])
    total_el = np.sum(np.abs(el_steps))
    if "short_stabilization" in indir:
        # Isolate the steps that last less than 12 minutes
        print("Assuming stabilization time is 5 minutes")
        in_observation = t_steps <= (12 * 60) / 86400
    else:
        # stabilization time is 30 minutes, calibration break is 5 minutes
        print("Assuming stabilization time is 30 minutes")
        in_observation = t_steps <= (37 * 60) / 86400
    obs_el = np.sum(np.abs(el_steps[in_observation]))
    steps = np.abs(el_steps[in_observation])
    large = steps > 1
    obs_el_large = np.sum(steps[large])
    average_throw = np.mean(np.abs(throws))

    average_rate = 1.0 / np.cos(np.radians(average_el))

    print("Schedule time:        {:.3f} days".format(schedule_time))
    print("Integration time:     {:.3f} days".format(integration_time))
    print("Observing efficiency: {:.3f} %".format(integration_time / schedule_time * 100))
    print("Integration time (Horizontal):     {:.3f} days".format(integration_time_horizontal))
    print("Observing efficiency (Horizontal): {:.3f} %".format(integration_time_horizontal / schedule_time * 100))
    print("Average elevation:    {:.3f} deg".format(average_el))
    print("Average throw:        {:.3f} deg".format(average_throw))
    print("Average Az-rate:      {:.3f} deg (for 1 deg/s on sky)".format(average_rate))
    print("Elevation travelled:  {:.3f} deg".format(total_el))
    print("Elevation travelled:  {:.3f} deg (during observation)".format(obs_el))
    print("Elevation travelled:  {:.3f} deg (requiring stabilization)".format(obs_el_large))

ax3.legend(loc="best")
fig.savefig("schedules.png")
plt.close()
