import numpy as np
import os
import sys
import matplotlib.pyplot as plt

import healpy as hp

# from toast_planck.utilities import plug_holes

plt.style.use("classic")

dec_min = -70
dec_max = 70
dec_step = 20
dec_overlap = 5
ra_step = 10
ra_overlap = 5

nbin = 10
# fwhm, tlim, plim = 10, 300 * 1e-6, 1 * 1e-7
fwhm, tlim, plim = 1, 500 * 1e-6, 2 * 1e-7
lmax = 512

wpatch = 10
hpatch = 20

south_tier_ranges = [(60, 30), (90, 10), (100, -65)]  # Tier 1, 2, 3
north_tier_ranges = [(165, 135), (175, 115), (230, 110)]  # Tier 1, 2, 3


def plug_holes(m, verbose=False, in_place=True, nest=False):
    """
    Use simple downgrading to derive estimates of the missing pixel values
    """
    nbad_start = np.sum(np.isclose(m, hp.UNSEEN))

    if nbad_start == m.size:
        if verbose:
            print("plug_holes: All map pixels are empty. Cannot plug holes", flush=True)
        return

    if nbad_start == 0:
        return

    nside = hp.get_nside(m)
    npix = m.size
    if nest:
        mnest = m.copy()
    else:
        mnest = hp.reorder(m, r2n=True)

    lowres = mnest
    nside_lowres = nside
    bad = np.isclose(mnest, hp.UNSEEN)
    while np.any(bad) and nside_lowres > 1:
        nside_lowres //= 2
        lowres = hp.ud_grade(lowres, nside_lowres, order_in="NESTED")
        hires = hp.ud_grade(lowres, nside, order_in="NESTED")
        bad = np.isclose(mnest, hp.UNSEEN)
        mnest[bad] = hires[bad]

    nbad_end = np.sum(bad)

    if nbad_end != 0:
        mn = np.mean(mnest[np.logical_not(bad)])
        mnest[bad] = mn

    if not in_place:
        m = m.copy()
    if nest:
        m[:] = mnest
    else:
        m[:] = hp.reorder(mnest, n2r=True)

    if verbose and nbad_start != 0:
        print(
            "plug_holes: Filled {} missing pixels ({:.2f}%), lowest "
            "resolution was Nside={}.".format(
                nbad_start, (100.0 * nbad_start) // npix, nside_lowres
            )
        )
    return m


def get_tier(left, right):
    tier = 4
    for ranges in [south_tier_ranges, north_tier_ranges]:
        for i, lims in enumerate(ranges):
            if left <= lims[0] and right >= lims[1]:
                return i + 1
    return tier


south_left = 100
south_top = [
    # Tier 3
    -15,  # 100
    -15,
    # Tier 2
    -20,  # 90
    -20,
    -25,  # 80
    # Tier 1
    -25,
    -30,  # 70
    -30,
    -30,  # 60
    # Tier 2
    -35,
    -35,  # 50
    -35,
    -30,  # 40
    -35,
    -35,  # 30
    -35,
    -35,  # 20
    # Tier 3
    -30,
    -30,  # 10
    -25,
    -25,  # 0,
    -25,
    -20,  # -10
    -20,
    -15,  # -20
    -15,
    -15,  # -30
    -10,
    -5,  # -40
    0,
    5,  # -50
    5,
]

north_left = 230
north_top = [
    # Tier 3
    -20,  # 230
    -25,
    -25,  # 220
    -25,
    -25,  # 210
    -25,
    -25,  # 200
    -20,
    -10,  # 190
    0,
    10,  # 180
    10,
    # Tier 2
    10,  # 170
    # Tier 1
    15,
    10,  # 160
    5,
    10,  # 150
    10,
    # Tier 2
    15,  # 140
    15,
    20,  # 130
    20,
    # Tier 3
    20,  # 120
]
"""
# left, top, right, bottom
primary_south = [80, -30, 30, -50]
primary_north = [175, 10, 135, -10]

secondary_south = [95, -10, 15, -60]
secondary_north = [180, 40, 125, -20]

tertiary_south = [100, -20, -60, -40]
tertiary_north = [210, 20, 110, 0]
"""

# target_top = -30
# target_bottom = -50
# target_left = 60
# target_right = 330

nside = 512
npix = 12 * nside ** 2

# fn_sync = "/Users/reijo/data/dx12/dx12_030_map.fits"
# fn_dust = "/Users/reijo/data/dx12/dx12_353_map.fits"
fn_sync = "/Users/reijo/data/hfi_pipe/commander_2018/sky_model_070GHz_nside1024_quickpol_cfreq_zodi.fits"
fn_dust = "/Users/reijo/data/hfi_pipe/commander_2018/sky_model_143GHz_nside2048_quickpol_cfreq_zodi.fits"
fn_mask = "/Users/reijo/data/hfi_pipe/psmask_030_nest.fits"
fn_cmb = "/Users/reijo/data/hfi_pipe/COM_CMB_IQU-commander_1024_R2.02_full.fits"
mask = hp.read_map(fn_mask, nest=True)
nside_mask = hp.get_nside(mask)
cmb = hp.read_map(fn_cmb, nest=True)
nside_cmb = hp.get_nside(cmb)

plt.figure(figsize=[18, 8])
plt.suptitle("fwhm = {} deg, nbin = {}".format(fwhm, nbin))


def get_p(fn):
    m_gal = hp.read_map(fn, range(3), nest=True)
    m_gal = hp.ud_grade(m_gal, nside_cmb, order_in="NEST")
    m_gal[0] -= cmb
    m_gal = hp.ud_grade(m_gal, nside_mask, order_in="NEST")
    m_gal[:, mask == 0] = hp.UNSEEN
    plug_holes(m_gal, nest=True, in_place=True, verbose=True)
    m_gal = hp.ud_grade(m_gal, 512, order_in="NEST", order_out="RING")
    m_gal = hp.smoothing(m_gal, fwhm=np.radians(fwhm), lmax=lmax, iter=0)
    i_gal = m_gal[0] - np.amin(m_gal[0])
    p_gal = np.sqrt(m_gal[1] ** 2 + m_gal[2] ** 2)
    i_gal[i_gal < tlim] = 0
    p_gal[p_gal < plim] = 0
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


isync_gal, psync_gal, isync_equ, psync_equ = get_p(fn_sync)
idust_gal, pdust_gal, idust_equ, pdust_equ = get_p(fn_dust)


def plot_p(maps, sub, title, nbin, coord, gr=30):
    sorted_maps = [np.sort(m) for m in maps]
    lims = np.arange(nbin) / nbin
    p = np.zeros(npix)
    for ilim, lim in enumerate(lims):
        for m, msorted in zip(maps, sorted_maps):
            val = msorted[int(npix * lim)]
            p[m > val] = lim
    hp.mollview(p, xsize=2400, sub=sub, title=title, coord=coord)
    hp.graticule(gr)
    return p


plot_p([isync_gal], [2, 4, 1], "100GHz I", nbin, "G")
plot_p([psync_gal], [2, 4, 2], "100GHz P", nbin, "G")
plot_p([idust_gal], [2, 4, 3], "150GHz I", nbin, "G")
plot_p([pdust_gal], [2, 4, 4], "150GHz P", nbin, "G")

plot_p([isync_equ], [2, 4, 5], "100GHz I", nbin, "C")
plot_p([psync_equ], [2, 4, 6], "100GHz P", nbin, "C")
plot_p([idust_equ], [2, 4, 7], "150GHz I", nbin, "C")
plot_p([pdust_equ], [2, 4, 8], "150GHz P", nbin, "C")

plt.savefig("priority_by_pixel_fwhm{:02}_nbin{:02}.comp.png".format(fwhm, nbin))

plt.figure(figsize=[18, 6])
plt.suptitle("fwhm = {} deg, nbin = {}".format(fwhm, nbin))

plot_p(
    [isync_gal, psync_gal, idust_gal, pdust_gal],
    [1, 2, 1],
    "100+150GHz",
    nbin,
    "G",
    gr=15,
)
pmap = plot_p(
    [isync_equ, psync_equ, idust_equ, pdust_equ],
    [1, 2, 2],
    "100+150GHz",
    nbin,
    "C",
    gr=15,
)


def plot_patch(patch, tier):
    left, top, right, bottom = patch
    color, lw = [("red", 3), ("black", 2), ("grey", 2)][tier - 1]
    lon = [left, left, right, right, left]
    lat = [bottom, top, top, bottom, bottom]
    n = len(lon)
    x = np.arange(n) / (n - 1)
    xfull = np.linspace(0, 1, 100 * n)
    lonfull = np.interp(xfull, x, lon)
    latfull = np.interp(xfull, x, lat)
    hp.projplot(
        lonfull,
        latfull,
        "-",
        threshold=1,
        lonlat=True,
        color=color,
        lw=lw,
        alpha=0.8,
        coord="C",
        zorder=(4 - tier) * 100,
    )


def add_patch(fout, patch, tier, mult=1):
    ra_stop, dec_start, ra_start, dec_stop = patch
    if ra_start < 0:
        ra_start += 360
    if ra_stop < 0:
        ra_stop += 360
    name = "Tier{}DEC{:+04}..{:+04}_RA{:+04}..{:+04}".format(
        tier, dec_start, dec_stop, ra_start, ra_stop
    )
    priority = [1, 1000, 1000000][tier - 1] * mult
    fout.write("--patch\n")
    fout.write(
        "{},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f}\n".format(
            name, priority, ra_start, dec_stop, ra_stop, dec_start
        )
    )


with open("patches.txt", "w") as fout:
    left = south_left
    for top in south_top:
        tier = get_tier(left, left - wpatch)
        patch = [left, top, left - wpatch, top - hpatch]
        plot_patch(patch, tier)
        add_patch(fout, patch, tier, mult=1)
        left -= 5

    left = north_left
    for top in north_top:
        tier = get_tier(left, left - wpatch)
        patch = [left, top, left - wpatch, top - hpatch]
        plot_patch(patch, tier)
        add_patch(fout, patch, tier, mult=10)
        left -= 5

plt.savefig("priority_by_pixel_fwhm{:02}_nbin{:02}.full.png".format(fwhm, nbin))
