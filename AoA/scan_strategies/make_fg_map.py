import numpy as np
import os
import sys
import matplotlib.pyplot as plt

import healpy as hp

from plug_holes import plug_holes

# Number of discrete colors to use in plotting
# intensity and polarization amplitude

nbin = 10
# fwhm, tlim, plim = 10, 300 * 1e-6, 1 * 1e-7
fwhm, tlim, plim = 1, 500 * 1e-6, 2 * 1e-7
lmax = 512

nside = 512
npix = 12 * nside ** 2

hfi_pipe = "/global/cfs/cdirs/cmb/data/planck2020/npipe/aux"
sky_model_cache = "/global/cfs/cdirs/cmb/data/planck2020/npipe/npipe6v20_sim/skymodel_cache"
fn_sync = f"{sky_model_cache}/sky_model_070GHz_nside1024_quickpol_cfreq_zodi.fits"
fn_dust = f"{sky_model_cache}/sky_model_143GHz_nside2048_quickpol_cfreq_zodi.fits"
fn_mask = f"{hfi_pipe}/psmask_030_nest.fits"
fn_cmb = f"{hfi_pipe}/COM_CMB_IQU-commander_1024_R2.02_full.fits"
mask = hp.read_map(fn_mask, nest=True)
nside_mask = hp.get_nside(mask)
cmb = hp.read_map(fn_cmb, nest=True)

plt.figure(figsize=[18, 8])
plt.suptitle("fwhm = {} deg, nbin = {}".format(fwhm, nbin))


def get_p(fn, cmb, psmask, fwhm, tlim, plim):
    """
    Return the intensity and polarization amplitude
    of the file in Galactic and Equatorial coordinates

    The map is first CMB-subtracted.  Outputs are
    smoothed to `fwhm` [degrees]
    """
    nside_cmb = hp.get_nside(cmb)
    m_gal = hp.read_map(fn, range(3), nest=True)
    # Subtract CMB
    m_gal = hp.ud_grade(m_gal, nside_cmb, order_in="NEST")
    m_gal[0] -= cmb
    # Mask out point sources
    nside_mask = hp.get_nside(psmask)
    m_gal = hp.ud_grade(m_gal, nside_mask, order_in="NEST")
    m_gal[:, mask == 0] = hp.UNSEEN
    plug_holes(m_gal, nest=True, in_place=True, verbose=True)
    # Smooth to FWHM
    m_gal = hp.ud_grade(m_gal, 512, order_in="NEST", order_out="RING")
    m_gal = hp.smoothing(m_gal, fwhm=np.radians(fwhm), lmax=lmax, iter=0)
    # Set zero point
    i_gal = m_gal[0] - np.amin(m_gal[0])
    # Measure polarization amplitude
    p_gal = np.sqrt(m_gal[1] ** 2 + m_gal[2] ** 2)
    # Saturate the color scale
    i_gal[i_gal < tlim] = 0
    p_gal[p_gal < plim] = 0
    # Rotate and interpolate to equatorial coordinates
    i_equ = np.zeros(npix)
    p_equ = np.zeros(npix)
    rot = hp.Rotator(coord="CG")
    buflen = 1024
    for pix_start in range(0, npix, buflen):
        pix_stop = min(npix, pix_start + buflen)
        pix = np.arange(pix_start, pix_stop)
        vec_equ = hp.pix2vec(nside, pix)
        vec_gal = rot(vec_equ)
        theta, phi = hp.vec2dir(vec_gal)
        ival = hp.get_interp_val(i_gal, theta, phi)
        pval = hp.get_interp_val(p_gal, theta, phi)
        i_equ[pix] = ival
        p_equ[pix] = pval
    return i_gal, p_gal, i_equ, p_equ


isync_gal, psync_gal, isync_equ, psync_equ = get_p(fn_sync, cmb, mask, fwhm, tlim, plim)
idust_gal, pdust_gal, idust_equ, pdust_equ = get_p(fn_dust, cmb, mask, fwhm, tlim, plim)


def plot_p(maps, sub, title, nbin, coord, gr=30):
    """
    Plot the given map using only `nbin` distinct colors
    If `maps` contains several maps, the plotting color is
    chosen based on the most intense map: if any of the given
    """
    sorted_maps = [np.sort(m) for m in maps]
    lims = np.arange(nbin) / nbin
    p = np.zeros(npix)
    for ilim, lim in enumerate(lims):
        for m, msorted in zip(maps, sorted_maps):
            val = msorted[int(npix * lim)]
            p[m > val] = lim
    hp.mollview(p, xsize=2400, sub=sub, title=title, coord=coord, cmap="bwr")
    hp.graticule(gr)
    return p


plot_p([isync_gal], [2, 4, 1], "070GHz I", nbin, "G")
plot_p([psync_gal], [2, 4, 2], "070GHz P", nbin, "G")
plot_p([idust_gal], [2, 4, 3], "150GHz I", nbin, "G")
plot_p([pdust_gal], [2, 4, 4], "150GHz P", nbin, "G")

plot_p([isync_equ], [2, 4, 5], "070GHz I", nbin, "C")
plot_p([psync_equ], [2, 4, 6], "070GHz P", nbin, "C")
plot_p([idust_equ], [2, 4, 7], "150GHz I", nbin, "C")
plot_p([pdust_equ], [2, 4, 8], "150GHz P", nbin, "C")

plt.savefig("priority_by_pixel_fwhm{:02}_nbin{:02}.comp.png".format(fwhm, nbin))

plt.figure(figsize=[18, 6])
plt.suptitle("fwhm = {} deg, nbin = {}".format(fwhm, nbin))

map_gal = plot_p(
    [isync_gal, psync_gal, idust_gal, pdust_gal],
    [1, 2, 1],
    "070+150GHz",
    nbin,
    "G",
    gr=15,
)

map_equ = plot_p(
    [isync_equ, psync_equ, idust_equ, pdust_equ],
    [1, 2, 2],
    "070+150GHz",
    nbin,
    "C",
    gr=15,
)

hp.write_map("fg_map_gal.fits", map_gal, nest=False, coord="G")
hp.write_map("fg_map_equ.fits", map_equ, nest=False, coord="C")

plt.savefig("priority_by_pixel_fwhm{:02}_nbin{:02}.full.png".format(fwhm, nbin))
