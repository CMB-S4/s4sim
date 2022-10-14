from glob import glob
import os
import sys

import matplotlib.pyplot as plt
import numpy as np


do_plot = False

if do_plot:
    fig = plt.figure(figsize=[18, 12])
    nrow, ncol = 2, 2
    ax1 = fig.add_subplot(nrow, ncol, 1)
    ax2 = fig.add_subplot(nrow, ncol, 2)
    ax3 = fig.add_subplot(nrow, ncol, 3)

    colors = ["tab:blue", "tab:orange", "tab:green", "tab:pink"]

if len(sys.argv) == 1:
    print("Usage: {sys.argv[0]} <input directory> [schedule time days]")
indir = sys.argv[1]
if len(sys.argv) == 3:
    schedule_days = float(sys.argv[2])
else:
    schedule_days = None
fnames = sorted(glob(indir + "/*txt"))

for i, fname in enumerate(fnames):
    if do_plot:
        color = colors[i]
    print("\n{}".format(fname))

    with open(fname, "r") as fin:
        fin.readline()
        header = fin.readline()
    try:
        site, telescope, lat, lon, alt = header.split()
    except:
        print(f"ERROR: failed to parse {fname}", flush=True)
        continue
    if np.abs(float(lat)) > 85:
        pole_site = True
    else:
        pole_site = False

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

    if lengths.size == 0:
        print(f"No scans found in {fname}")
        continue

    if lengths.size == 1:
        print(f"Only one scan in {fname}")
        continue

    last_el = None
    line = 3
    for start, stop, az_min, az_max, el, name in zip(
        arr["starts"], arr["stops"], arr["az_mins"], arr["az_maxs"], arr["elevations"], arr["names"],
    ):
        line += 1
        if not pole_site and az_min < 180 and az_max > 180:
            print(
                f"WARNING: Found a mixed Rising/Setting scan: {az_min} - {az_max}. line = {line}",
                flush=True,
            )
        if last_el is None:
            label = os.path.basename(fname).replace(".txt", "")
        else:
            label = None
        if do_plot:
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

    if schedule_days is not None:
        print(f"Overriding schedule_time = {schedule_time} with {schedule_days}")
        schedule_time = schedule_days

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
    print("Elevation travelled:  {:.3f} deg (requiring stabilization)".format(obs_el_large), flush=True)

if do_plot:
    ax3.legend(loc="best")
    fig.savefig("schedules.png")
    plt.close()
