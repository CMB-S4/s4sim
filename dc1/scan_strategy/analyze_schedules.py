import os
import sys

import numpy as np


for fname in [
        "pole_sat/schedules/pole_schedule_sat.txt",
        "chile_lat/schedules/chile_schedule_lat.txt", 
        "pole_lat/schedules/pole_schedule_lat.txt",
        ]:
    arr = np.genfromtxt(fname, skip_header=3, unpack=True)
    starts, stops = arr[4], arr[5]
    lens = stops - starts
    
    gaps = starts[1:] - stops[:-1]
    observing_breaks = np.argwhere(gaps > 1.0).ravel()
    if len(observing_breaks) == 0:
        nbreak = 0
        break_len = 0
    else:
        nbreak = len(observing_breaks)
        break_len = np.sum(gaps[observing_greaks])
    
    good = lens < (10 / 24)
    total_len = stops[-1] - starts[0]
    good_len = np.sum(lens[good])
    bad_len = np.sum(lens[np.logical_not(good)])

    print(
        f"{fname:>45} : total schedule = {total_len:.2f} days, "
        f"{nbreak} observing_breaks = {break_len:.2f} days, "
        f"science = {good_len:.2f} days ({100 * good_len / total_len:.2f}%), "
        f"calibration = {bad_len:.2f} days ({100 * bad_len / total_len:.2f}%)"
    )
