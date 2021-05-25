import os
import sys

import healpy as hp
import numpy as np

import s4sim.hardware as hardware


telescope = "LAT0"
ntele = 2

ndet_ref = 5832
nday_ref = 10
band_ref = "MFL1"

nday = 7 * 365

print("Simulating hardware", flush=True)
hw = hardware.get_example()

print("Simulating detectors", flush=True)
hw.data["detectors"] = hardware.sim_telescope_detectors(hw, telescope)
bdata_ref = hw.data["bands"][band_ref]

for el in 30, 35, 40, 45, 50, 55, 60:
    sin_el = np.sin(np.radians(el))
    net_ref = bdata_ref["NET"] * (bdata_ref["A"] / sin_el + bdata_ref["C"])

    fname_in = f"out-experimental/00000000/chile_noise_LAT_{band_ref}_el{el}_filtered_telescope_all_time_all_wcov.fits"
    print(f"Reading {fname_in}", flush=True)
    w = hp.read_map(fname_in, None, verbose=False)
    bad = np.logical_or(w[0] == hp.UNSEEN, w[0] == 0)
    # Discard worst pixels. This limit is arbitrary.
    not_bad = np.logical_not(bad)
    good = not_bad.copy()
    wsorted = np.sort(w[0, good])
    lim = wsorted[int(0.99 * wsorted.size)]
    not_bad[w[0] > lim] = False
    
    lim = 4 * np.median(w[0, good])
    good[w[0] > lim] = False
    fsky = np.sum(good) / good.size
    fsky_raw = 1 - np.sum(bad) / good.size
    nside = hp.get_nside(w)
    pix_area = hp.nside2pixarea(nside, degrees=True)
    pix_side = np.sqrt(pix_area) * 60  # in arc min

    print(f"el = {el}, fsky(99%) = {fsky_raw:8.3} (fsky(good) = {fsky:8.3})")

    for band in "LFL1", "LFL2", "MFL1", "MFL2", "HFL1", "HFL2":
        # Observing efficiency from
        #   https://docs.google.com/spreadsheets/d/1jR9gSsJ0w1dEO5Jb_URlD3SWYtgFtwBgB3W88p6puo0/edit?usp=sharing
        obs_eff = {
            "LFL1" : 0.25,
            "LFL2" : 0.25,
            "MFL1" : 0.25,
            "MFL2" : 0.25,
            "HFL1" : 0.22,
            "HFL2" : 0.22,
        }[band]

        #print(f"Selecting band = {band}", flush=True)
        hw_band = hw.select(tubes=None, match={"band" : band})
        ndet = len(hw_band.data["detectors"])

        bdata = hw_band.data["bands"][band]
        net =  bdata["NET"] * (bdata["A"] / sin_el + bdata["C"])

        w_band = w * ndet_ref / ndet * nday_ref / nday * net / net_ref / ntele / obs_eff
        w_band[:, bad] = 0

        depth = np.vstack([w_band[0], w_band[3], w_band[5]]) ** .5 * 1e6 * pix_side
        depth[:, bad] = hp.UNSEEN

        fname_out = f"out-experimental/00000000/noise_depth_LAT_{band}_el{el}.fits"
        #if not os.path.isfile(fname_out):
        hp.write_map(fname_out, depth, overwrite=True)

        depth_I_raw = np.mean(depth[0, not_bad])
        depth_Q_raw = np.mean(depth[1, not_bad])
        depth_U_raw = np.mean(depth[2, not_bad])
        depth_I = np.mean(depth[0, good])
        depth_Q = np.mean(depth[1, good])
        depth_U = np.mean(depth[2, good])
        print(
            f"  band = {band}, ndet = {ndet:8}, "
            f"IQU depth(99%) = {depth_I_raw:10.3f}, {depth_Q_raw:10.3f}, {depth_U_raw:10.3f}, ",
            f"IQU depth(good) = {depth_I:10.3f}, {depth_Q:10.3f}, {depth_U:10.3f}",
            flush=True,
        )
