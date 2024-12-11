import dateutil
import os
import sys

import ephem
import matplotlib.pyplot as plt
import numpy as np

import toast.coordinates
import toast.utils
import toast.weather


np.random.seed(673653982)
realization = 0

pwv_limits = [2, 3]
fname_in = "schedule_sat.sun90.txt"
# fname_in = "schedule_sat.sun45.txt"

first_time = dateutil.parser.parse(f"2030-01-01 00:00:00 00:00:00 +0000")
break_start = dateutil.parser.parse(f"2030-01-01 00:00:00 00:00:00 +0000")
break_stop = dateutil.parser.parse(f"2030-04-01 00:00:00 00:00:00 +0000")
break_frac = 19 / 89  # adjust pwv_limit during the break to meet this efficiency

for pwv_limit in pwv_limits:
    print(f"{fname_in}, PWV limit = {pwv_limit}")
    fname_out = fname_in.replace(".txt", f".{pwv_limit}mm.txt")

    with open(fname_in, "r") as schedule_in:
        header = []
        pwvs = []
        times = []
        lengths = []
        sun_els = []
        lines = []
        break_flags = []
        break_pwvs = []
        break_times = []
        break_lengths = []
        break_sun_els = []
        break_lines = []
        read_header = True
        weather = None
        observer = ephem.Observer()
        sun = ephem.Sun()
        for line in schedule_in:
            if line.startswith("#"):
                header.append(line)
                continue
            if read_header:
                site_name, telescope_name, site_lat, site_lon, site_alt = line.split()
                observer.lon = site_lon
                observer.lat = site_lat
                observer.elevation = float(site_alt)
                observer.epoch = "2000"
                observer.temp = 0  # in Celcius
                observer.compute_pressure()
                site_uid = toast.utils.name_UID(site_name)
                header.append(line)
                read_header = False
                continue

            parts = line.split()
            name = parts[7]
            if name == "CALIBRATION_BREAK":
                continue

            start_date, start_time, stop_date, stop_time = parts[:4]
            start_time = dateutil.parser.parse(f"{start_date} {start_time} +0000")
            stop_time = dateutil.parser.parse(f"{stop_date} {stop_time} +0000")
            mid_time = start_time + 0.5 * (stop_time - start_time)
            if weather is None:
                weather = toast.weather.SimWeather(
                    time=mid_time,
                    name="atacama",
                    site_uid=site_uid,
                    realization=realization,
                )
            else:
                weather.set(time=mid_time, realization=realization, site_uid=site_uid)
            # Get Sun's elevation at mid_time
            observer.date = toast.coordinates.to_DJD(mid_time.timestamp())
            sun.compute(observer)
            sun_el = np.degrees(sun.alt)

            during_break = break_start < start_time and start_time < break_stop
            if during_break:
                break_pwvs.append(weather.pwv.to_value("mm"))
                break_times.append((mid_time - first_time).total_seconds())
                break_lengths.append((stop_time - start_time).total_seconds())
                break_sun_els.append(sun_el)
                break_lines.append(line)
            else:
                pwvs.append(weather.pwv.to_value("mm"))
                times.append((mid_time - first_time).total_seconds())
                lengths.append((stop_time - start_time).total_seconds())
                sun_els.append(sun_el)
                lines.append(line)

    pwvs = np.array(pwvs)
    times = np.array(times)
    lengths = np.array(lengths)
    sun_els = np.array(sun_els)
    break_flags = np.array(break_flags)
    break_pwvs = np.array(break_pwvs)
    break_times = np.array(break_times)
    break_lengths = np.array(break_lengths)
    break_sun_els = np.array(break_sun_els)

    ind = np.argsort(pwvs)
    sorted_pwvs = pwvs[ind]
    sorted_lengths = lengths[ind]

    n = np.searchsorted(sorted_pwvs, pwv_limit)
    fyear = np.sum(lengths[:n]) / np.sum(lengths)
    print(f"f_year for PWV limit = {pwv_limit} mm is {fyear}")

    ind = np.argsort(break_pwvs)
    sorted_break_pwvs = break_pwvs[ind]
    sorted_break_lengths = break_lengths[ind]
    sorted_break_sun_els = break_sun_els[ind]
    # First cut observations with high PWV.
    # If the break_frac is not met, start cutting highest PWV
    # daytime observations
    target_length = np.sum(sorted_break_lengths) * break_frac
    good = sorted_break_pwvs < pwv_limit
    n = len(sorted_break_lengths)
    daytime_pwv_limit = pwv_limit
    for i in range(n - 1, -1, -1):
        j = n - 1 - i
        if np.sum(sorted_break_lengths[good]) < target_length:
            break
        if good[i] and sorted_break_sun_els[i] > 0:
            daytime_pwv_limit = sorted_break_pwvs[i]
            good[i] = False
    pwv_limit_break = daytime_pwv_limit
    frac_eff = np.sum(sorted_break_lengths[good]) / np.sum(sorted_break_lengths)
    print(
        f"To meet eff = {break_frac:.3f}, break requires "
        f"PWV limit = {pwv_limit_break:.3f}. "
        f"eff = {frac_eff:.3f}"
    )

    # Plot

    fig = plt.figure(figsize=[18, 12])
    nrow, ncol = 3, 1

    ax = fig.add_subplot(nrow, ncol, 1)
    ax.plot(times / 86400, pwvs, '.', label=f"{pwvs.size}")
    ax.plot(break_times / 86400, break_pwvs, '.', label=f"{break_pwvs.size}")
    ax.axhline(pwv_limit, color='k', linestyle="--")
    ax.set_xlabel("DOY")
    ax.set_ylabel("PWV [mm]")
    ax.set_title(f"All {len(lengths) + len(break_lengths)} observations")
    ax.set_xlim([0, 366])
    plt.legend(loc="best")

    ax = fig.add_subplot(nrow, ncol, 2)
    good = pwvs < pwv_limit
    good_break = break_pwvs < pwv_limit
    good_break[np.logical_and(
        break_sun_els > 0, break_pwvs > daytime_pwv_limit
    )] = False
    ngood = np.sum(good)
    ngood_break = np.sum(good_break)
    ax.plot(
        times[good] / 86400, pwvs[good],
        '.',
        label=f"{ngood} : PWV < {pwv_limit}",
    )
    ax.plot(
        break_times[good_break] / 86400,
        break_pwvs[good_break],
        '.',
        label=f"{ngood_break} : PWV < {pwv_limit}/{daytime_pwv_limit:.3f}",
    )
    ax.set_xlabel("DOY")
    ax.set_ylabel("PWV [mm]")
    ax.set_title(f"observations after PWV cuts")
    ax.set_xlim([0, 366])
    plt.legend(loc="best")

    # Compute observing efficiencies

    all_times = np.hstack([break_times, times])
    all_lengths = np.hstack([break_lengths, lengths])
    all_sun_els = np.hstack([break_sun_els, sun_els])
    good_times = np.hstack([break_times[good_break], times[good]])
    good_lengths = np.hstack([break_lengths[good_break], lengths[good]])
    t_month = 30.5 * 86400
    x, y0, y = [], [], []  # for plotting
    print("Monthly observing efficiencies after PWV cut.")
    for month in range(12):
        tstart = month * t_month
        tstop = tstart + t_month
        ind = np.logical_and(good_times > tstart, good_times < tstop)
        frac = np.sum(good_lengths[ind]) / t_month
        # Uncut, for plotting
        ind0 = np.logical_and(all_times > tstart, all_times < tstop)
        frac0 = np.sum(all_lengths[ind0]) / t_month
        x += [tstart / 86400, tstop / 86400]
        y0 += [frac0, frac0]
        y += [frac, frac]
        print(
            f"{month + 1:02} : Total: {frac:.3f} "
            f"(ideal = {frac0:.3f}, PWV cut = {frac / frac0:.3f})")

    ax = fig.add_subplot(nrow, ncol, 3)
    ax.plot(x, y, label="w/ PWV limit")
    ax.plot(x, y0, label="w/o PWV limit")
    ax.plot(x, np.array(y) / np.array(y0), label="PWV cut")
    ax.set_ylim([0, 1.1])
    ax.axhline(1, color="k")
    ax.legend(loc="best")
    ax.set_xlabel("DOY")
    ax.set_ylabel("Obs. efficiency")
    ax.set_xlim([0, 366])

    # Save the plot

    fig.subplots_adjust(hspace=0.4)
    fname_plot = fname_in.replace(".txt", ".png")
    if fname_plot == schedule_in:
        fname_plot += ".png"
    fname_plot = fname_plot.replace(".png", f".{pwv_limit}mm.png")
    fig.savefig(fname_plot)
    print(f"Wrote {fname_plot}")

    # write out the observations that survive the cuts

    obs_time_day = 0
    obs_time_night = 0
    with open(fname_out.replace(".txt", ".season.txt"), "w") as schedule_out:
        for line in header:
            schedule_out.write(line)
        for pwv, length, line, flag, sun_el in zip(pwvs, lengths, lines, good, sun_els):
            if not flag:
                continue
            if sun_el > 0:
                obs_time_day += length
            else:
                obs_time_night += length
            schedule_out.write(line)

    obs_time_break_day = 0
    obs_time_break_night = 0
    with open(fname_out.replace(".txt", ".break.txt"), "w") as schedule_out:
        for line in header:
            schedule_out.write(line)
        for pwv, length, line, flag, sun_el in zip(
                break_pwvs,
                break_lengths,
                break_lines,
                good_break,
                break_sun_els,
        ):
            if not flag:
                continue
            if sun_el > 0:
                obs_time_break_day += length
            else:
                obs_time_break_night += length
            schedule_out.write(line)

    break_length = (break_stop - break_start).total_seconds()
    year_length = 365 * 86400

    year_length_day = year_length * 0.5
    year_length_night = year_length * 0.5
    f_field_total = np.sum(all_lengths) / year_length
    f_field_total_day = np.sum(all_lengths[all_sun_els > 0]) / year_length_day
    f_field_total_night = np.sum(all_lengths[all_sun_els < 0]) / year_length_night
    
    f_field_season = np.sum(lengths) / (year_length - break_length)
    f_field_break = np.sum(break_lengths) / break_length

    print(f"f_field(total)  = {f_field_total:.3f} (day: {f_field_total_day:.3f}, night: {f_field_total_night:.3f})")
    print(f"f_field(season) = {f_field_season:.3f}")
    print(f"f_field(break)  = {f_field_break:.3f}")

    day = 86400
    print(
        f"Surviving observing time (season): "
        f"{(obs_time_day + obs_time_night) / day:.3f} days "
        f"(day: {obs_time_day / 86400:.3f}, night: {obs_time_night / day:.3f})"
    )
    print(
        f"Surviving observing time (break): "
        f"{(obs_time_break_day + obs_time_break_night) / day:.3f} days "
        f"(day: {obs_time_break_day / 86400:.3f}, night: {obs_time_break_night / day:.3f})"
    )
