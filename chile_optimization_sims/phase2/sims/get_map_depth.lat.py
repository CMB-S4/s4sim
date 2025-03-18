import glob
import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np


bands = [20, 30, 40, 90, 150, 220, 280]
add_45 = False
add_wide = True
just_wide = False

def fskies(hit):
    good = hit > 0
    dhit = hit[good].astype(float)
    moment0 = np.sum(good) / npix
    moment1 = np.sum(dhit) / npix
    moment2 = np.sum(dhit ** 2) / npix
    moment4 = np.sum(dhit ** 4) / npix
    fraw = moment0
    fnoise = moment1 ** 2 / moment2
    fsignal = moment2 ** 2 / moment4
    return fraw, fnoise, fsignal

def invert_map(m):
    minv = np.zeros_like(m)
    good = m != 0
    minv[good] = 1 / m[good]
    return minv


band = 90
satcov = hp.read_map(
    glob.glob(f"scaled_outputs/sun90max_f{band:03}_*years_cov.fits")[0]
)
if add_45:
    cov2 = hp.read_map(
        glob.glob(f"scaled_outputs/sun45_supplement_f{band:03}_*years_cov.fits")[0]
    )
    satcov = invert_map(invert_map(satcov) + invert_map(cov2))
satgood = satcov != 0
satinvcov = np.zeros_like(satcov)
satinvcov[satgood] = 1 / satcov[satgood]


# for nlat in 2, 3, 4, 5:
for nlat in 2,:
    if nlat > 3 and just_wide:
        break
    for band in bands:
        print(band)
        if nlat == 2:
            # Reduced configuration, assume 10 LAT years of delensing
            cov = hp.read_map(
                glob.glob(f"scaled_outputs/lat_delensing_max_f{band:03}_16years_cov.fits")[0]
            )
            cov *= 16 / 10
        else:
            cov = hp.read_map(
                glob.glob(f"scaled_outputs/lat_delensing_max_f{band:03}_{6+(nlat-2)*10}years_cov.fits")[0]
            )
        if add_wide or just_wide:
            cov2 = hp.read_map(
                glob.glob(f"scaled_outputs/lat_wide_f{band:03}_14years_cov.fits")[0]
            )
            if just_wide:
                cov = cov2
            else:
                cov = invert_map(invert_map(cov) + invert_map(cov2))
        good = cov != 0
        invcov = np.zeros_like(cov)
        invcov[good] = 1 / cov[good]

        nside = hp.get_nside(cov)
        npix = cov.size
        lon, lat = hp.pix2ang(nside, np.arange(npix), lonlat=True)
        pix_area = hp.nside2pixarea(nside, degrees=True) * 3600  # arcmin^2


        nrow, ncol = 1, 3
        fig = plt.figure(figsize=[ncol * 6, nrow * 4])
        if add_wide:
            fig.suptitle(f"LAT {band}GHz delens + wide, nlat={nlat}")
        else:
            fig.suptitle(f"LAT {band}GHz delens, nlat={nlat}")

        iplot = 0
        for field in "all", "south", "north":
            iplot += 1
            hits = invcov.copy()
            if field == "all":
                mask = hp.read_map("cmbs4_phase2_mask_south.fits", dtype=bool)
                mask += hp.read_map("cmbs4_phase2_mask_north.fits", dtype=bool)
                if not just_wide:
                    hits *= mask
            elif field == "south":
                hits *= hp.read_map("cmbs4_phase2_mask_south.fits", dtype=bool)
            elif field == "north":
                hits *= hp.read_map("cmbs4_phase2_mask_north.fits", dtype=bool)
            frac = np.sum(hits) / np.sum(invcov)
            # hits[satinvcov == 0] = 0
            good_field = hits != 0
            rms = np.sqrt(cov) * good_field
            depth_I = rms * np.sqrt(pix_area) * 1e6  # uK.arcmin
            depth = np.sqrt(2) * depth_I
            rms_eff = np.sqrt(np.sum(invcov[good_field]) / np.sum(invcov[good_field]**2)) * np.sqrt(pix_area) * 1e6 * np.sqrt(2)  # uK.arcmin (pol)
            # rms_eff = np.sqrt(np.sum(cov[good_field] / satcov[good_field]**2) / np.sum(satcov[good_field]**2)) * np.sqrt(pix_area) * 1e6 * np.sqrt(2)  # uK.arcmin (pol)
            depth[depth == 0] = hp.UNSEEN
            fraw, fnoise, fsignal = fskies(hits)
            # fraw, fnoise, fsignal = fskies(satinvcov * good_field)
            vmin = np.amin(depth[good_field])
            vmax = 2 * vmin
            hp.mollview(
                depth,
                title=f"{field} : fsky={fraw:.3f}/{fnoise:.3f}/{fsignal:.3f} ; depth={rms_eff:.3f}, frac={frac:0.3f}",
                sub=[nrow, ncol, iplot],
                cmap="magma",
                unit=r"$\mu$K-arcmin",
                min=vmin,
                max=vmax,
            )

        if just_wide:
            fig.savefig(f"lat_depth_{band:03}.wide.png")
        elif add_wide:
            fig.savefig(f"lat_depth_{band:03}.delens+wide.nlat={nlat}.png")
        else:
            fig.savefig(f"lat_depth_{band:03}.delens.nlat={nlat}.png")
        plt.close()
