import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np


plot_bk = True

fname_pdust = "pdust.fits"
if os.path.isfile(fname_pdust):
    print(f"Loading {fname_pdust}")
    pdust = hp.read_map(fname_pdust)
else:
    fname_dust = "/global/cfs/cdirs/cmb/data/planck2020/npipe/npipe6v20/npipe6v20_353_map.fits"
    print(f"Loading {fname_dust}")
    dust = hp.read_map(fname_dust, None)
    print("Smoothing")
    sdust = hp.smoothing(dust, fwhm=np.radians(1), lmax=512, iter=0)
    pdust = np.sqrt(sdust[1]**2 + sdust[2]**2)
    print("Downgrading")
    pdust = hp.ud_grade(pdust, 512)
    print("Rotating")
    rot = hp.Rotator(coord="GC")
    pdust = rot.rotate_map_alms(pdust)
    print(f"Writing {fname_pdust}")
    hp.write_map(fname_pdust, pdust)

if plot_bk:
    fname_bk = "bk18_mask_largefield_cel_n0512.fits"
    print(f"Loading {fname_bk}")
    bkhits = hp.read_map(fname_bk)
    bkhits[np.isnan(bkhits)] = 0

for flavor in "sat_deep", "sat_deep_bk":
    fname = {
        "sat_deep" : "outputs/sat_deep/SAT_f090/mapmaker_invcov.fits",
        "sat_deep_bk" : "outputs/sat_deep_bk/SAT_f090/mapmaker_invcov.fits",
    }[flavor]
    print(f"Loading {fname}")
    invcov = hp.read_map(fname)
    # depth = np.sqrt(invcov)
    weight = invcov / np.amax(invcov)
    good = weight != 0
    weight[weight == 0] = hp.UNSEEN

    nrow, ncol = 1, 2
    fig = plt.figure(figsize=[6 * ncol, 4 * nrow])
    hp.mollview(weight, sub=[nrow, ncol, 1], cmap="magma", title="", cbar=False, min=0, max=1)

    def get_outline(weight, frac):
        good = np.logical_and(weight != 0, weight != hp.UNSEEN)
        sorted_weight = np.sort(weight[good])[::-1]
        cumulative_weight = np.cumsum(sorted_weight)
        tot = np.sum(weight[good])

        i = np.searchsorted(cumulative_weight, frac * tot)
        mask = (weight > sorted_weight[i]).astype(float)
        fsky = np.sum(mask) / mask.size
        print(f"{frac} of the survey is limited to fsky = {fsky:.3f}")
        fsky = np.sum(good) / good.size
        print(f"The survey is limited to fsky = {fsky:.3f}")

        smask = hp.smoothing(mask, fwhm=np.radians(2), lmax=256)
        p = 0.10
        smask[smask > 1 - p] = 0
        smask[smask > p] = 1
        smask[smask < p] = 0

        return smask

    outline10 = get_outline(weight, 0.25)
    outline50 = get_outline(weight, 0.50)
    outline_bk10 = get_outline(bkhits, 0.25)
    outline_bk50 = get_outline(bkhits, 0.50)

    hp.mollview(pdust, sub=[nrow, ncol, 2], cmap="magma", title="", cbar=False, max=1e-4)
    hp.mollview(outline10, sub=[nrow, ncol, 2], cmap="bwr", title="", cbar=False, reuse_axes=True, alpha=0.75 * outline10)
    hp.mollview(outline50, sub=[nrow, ncol, 2], cmap="bwr", title="", cbar=False, reuse_axes=True, alpha=0.50 * outline50)
    hp.mollview(outline_bk10, sub=[nrow, ncol, 2], cmap="plasma", title="", cbar=False, reuse_axes=True, alpha=0.75 * outline_bk10)
    hp.mollview(outline_bk50, sub=[nrow, ncol, 2], cmap="plasma", title="", cbar=False, reuse_axes=True, alpha=0.50 * outline_bk50)

    fname_plot = f"survey_weight_{flavor}.png"
    plt.savefig(fname_plot)
    print(f"Wrote {fname_plot}")

plt.show()
