import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np


#for tele in "chlat", "splat", "spsat":
for tele in "spsat",:
    rootdir = f"/global/cfs/cdirs/cmbs4/dc/dc0/mission/{tele}"

    # Plot all splits for one frequency and one component

    nsplit = 32
    nrow, ncol = 8, 4
    if tele == "spsat":
        band = "155"
    else:
        band = "150"
    map_type = "map03"  # T/P depth map

    fig = plt.figure(figsize=[ncol * 4, nrow * 3])

    iplot = 0
    for isplit in range(1, nsplit + 1):
        iplot += 1

        fname_in = os.path.join(
            rootdir,
            f"split{nsplit:02}",
            band,
            f"dc0_{tele}_t{nsplit:02}.{isplit:02}_{band}_{map_type}.fits",
        )
        print(f"Loading {fname_in}", flush=True)
        m = hp.read_map(fname_in)
        m *= 1e6
        m[m == 0] = hp.UNSEEN
        if tele == "chlat":
            vmin, vmax = 5, 15
        else:
            vmin, vmax = 2, 5
        hp.mollview(
            m,
            sub=[nrow, ncol, iplot],
            title=f"{isplit:2} / {nsplit}",
            cmap="magma",
            unit="$\mu$K$_\mathrm{CMB}$ arcmin",
            min=vmin,
            max=vmax,
            xsize=800,
        )

    fig.tight_layout()

    fname_out = f"depths_{tele}_{nsplit}.png"
    fig.savefig(fname_out)
    print(f"Wrote {fname_out}", flush=True)

    fname_out = f"depths_{tele}_{nsplit}.pdf"
    fig.savefig(fname_out)
    print(f"Wrote {fname_out}", flush=True)

    nsplit = 32
    nrow, ncol = 8, 4
    # band = "150"
    complement = "1000"  # noise + atmosphere
    map_type = "map02"  # filter&bin IQU map

    fig = plt.figure(figsize=[ncol * 4, nrow * 3])

    iplot = 0
    for isplit in range(1, nsplit + 1):
        iplot += 1

        fname_in = os.path.join(
            rootdir,
            f"split{nsplit:02}",
            band,
            f"dc0_{tele}_t{nsplit:02}.{isplit:02}_{band}_{map_type}_c{complement}.fits",
        )
        print(f"Loading {fname_in}", flush=True)
        m = hp.read_map(fname_in)
        good = m != hp.UNSEEN
        m[good] *= 1e6
        if tele == "chlat":
            amp = 100
        else:
            amp = 30
        hp.mollview(
            m,
            sub=[nrow, ncol, iplot],
            title=f"{isplit:2} / {nsplit}",
            cmap="bwr",
            unit="$\mu$K$_\mathrm{CMB}$",
            min=-amp,
            max=amp,
        )

    fig.tight_layout()

    fname_out = f"splits_{tele}_{nsplit}.png"
    fig.savefig(fname_out)
    print(f"Wrote {fname_out}", flush=True)

    fname_out = f"splits_{tele}_{nsplit}.pdf"
    fig.savefig(fname_out)
    print(f"Wrote {fname_out}", flush=True)

    # Plot all components and all frequencies

    nsplit = 1
    isplit = 1
    map_type = "map02"  # filter&bin IQU map

    if tele == "chlat":
        bands = ["025", "040", "090", "150", "230", "280"]
    elif tele == "splat":
        bands = ["020" , "025", "040", "090", "150", "230", "280"]
    elif tele == "spsat":
        bands = ["025", "040", "085", "095", "145", "155", "230", "280"]
    complements = ["0001", "0010", "0100", "1000"]

    nrow = len(bands)
    ncol = len(complements)
    fig = plt.figure(figsize=[ncol * 4, nrow * 3])

    iplot = 0
    for band in bands:
        for complement in complements:
            iplot += 1

            complement_name = {
                "1000" : "atmosphere + noise",
                "0100" : "foregrounds",
                "0010" : "CMB lensing",
                "0001" : "unlensed CMB",
            }[complement]

            fname_in = os.path.join(
                rootdir,
                f"split{nsplit:02}",
                band,
                f"dc0_{tele}_t{nsplit:02}.{isplit:02}_{band}_{map_type}_c{complement}.fits",
            )
            print(f"Loading {fname_in}", flush=True)
            m = hp.read_map(fname_in)
            good = m != hp.UNSEEN
            m[good] *= 1e6
            amp = 100
            hp.mollview(
                m,
                sub=[nrow, ncol, iplot],
                title=f"{band} {complement_name}",
                cmap="bwr",
                unit="$\mu$K$_\mathrm{CMB}$",
                min=-amp,
                max=amp,
            )

    fig.tight_layout()

    fname_out = f"components_{tele}_{nsplit}.png"
    fig.savefig(fname_out)
    print(f"Wrote {fname_out}", flush=True)

    fname_out = f"components_{tele}_{nsplit}.pdf"
    fig.savefig(fname_out)
    print(f"Wrote {fname_out}", flush=True)
