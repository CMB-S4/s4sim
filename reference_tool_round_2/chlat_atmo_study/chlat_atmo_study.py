import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
import pylab

nside = 4096
npix = 12 * nside ** 2
lmax = 2 * nside
lmax_tf = 2048


def get_mask(fname_hits):
    fname_mask = "mask_" + os.path.basename(fname_hits)
    if os.path.isfile(fname_mask):
        mask = hp.read_map(fname_mask)
    else:
        hits = hp.read_map(fname_hits)
        good = hits > 0
        ngood = np.sum(good)
        sorted_hits = np.sort(hits[good])

        hit_lim = sorted_hits[np.int(ngood * .01)]
        mask = hits > hit_lim

        pix = np.arange(hits.size)
        lon, lat = hp.pix2ang(nside, pix, lonlat=True)
        lat_min = np.amin(lat[mask])
        lat_max = np.amax(lat[mask])

        mask = np.zeros(npix)
        tol = 10.0  # degrees
        mask[np.logical_and(lat_min + tol < lat, lat < lat_max - tol)] = 1
        mask = hp.smoothing(mask, fwhm=np.radians(3), lmax=2048)
        hp.write_map(fname_mask, mask)
    return mask


def map2cl(m, mask):
    m[0] = hp.remove_dipole(m[0])
    m[m == hp.UNSEEN] = 0
    fsky = np.sum(mask) / mask.size
    cl = hp.anafast(m * mask, lmax=lmax, iter=0) / fsky
    return cl


def get_cl(fname_map, fname_hits=None):
    fname_cl = os.path.basename(fname_map).replace("bmap", "cl")
    if fname_cl == os.path.basename(fname_map):
        fname_cl = "cl_" + fname_cl
    if os.path.isfile(fname_cl):
        cl = hp.read_cl(fname_cl)
    else:
        mask = get_mask(fname_hits)
        m = hp.read_map(fname_map, None)
        cl = map2cl(m, mask)
        hp.write_cl(fname_cl, cl)
    return cl


rootdir = "/global/cscratch1/sd/keskital/s4sim/reference_tool_round_2/out/00000000/"
fname_hits = os.path.join(rootdir, "chile_noise_LAT_MFL2_filtered_telescope_all_time_all_hmap.fits")
fname_atmo00 = os.path.join(rootdir, "chile_atmosphere_2DorderNO_LAT_MFL2_filtered_telescope_all_time_all_bmap.fits")
fname_atmo0 = os.path.join(rootdir, "chile_atmosphere_LAT_MFL2_filtered_telescope_all_time_all_bmap.fits")
#fname_atmo1 = os.path.join(rootdir, "chile_atmosphere_2Xdet_LAT_MFL2_filtered_telescope_all_time_all_bmap.fits")
fname_atmo1 = os.path.join(rootdir, "chile_atmosphere_2Dorder2_LAT_MFL2_filtered_telescope_all_time_all_bmap.fits")
fname_atmo2 = os.path.join(rootdir, "chile_atmosphere_2Dorder5_LAT_MFL2_filtered_telescope_all_time_all_bmap.fits")
fname_atmo3 = os.path.join(rootdir, "chile_atmosphere_2Dorder10_LAT_MFL2_filtered_telescope_all_time_all_bmap.fits")
#fname_atmo4 = os.path.join(rootdir, "chile_atmosphere_2Dorder20_LAT_MFL2_filtered_telescope_all_time_all_bmap.fits")
fname_atmo4 = os.path.join(rootdir, "chile_atmosphere_2Dorder0_by_tube_LAT_MFL2_filtered_telescope_all_time_all_bmap.fits")
fname_atmo5 = os.path.join(rootdir, "chile_atmosphere_2Dorder1_by_tube_LAT_MFL2_filtered_telescope_all_time_all_bmap.fits")
fname_cmb00 = os.path.join(rootdir, "chile_cmb_unlensed_solardipole_2DorderNO_LAT_MFL2_filtered_telescope_all_time_all_bmap.fits")
fname_cmb0 = os.path.join(rootdir, "chile_cmb_unlensed_solardipole_LAT_MFL2_filtered_telescope_all_time_all_bmap.fits")
fname_cmb1 = os.path.join(rootdir, "chile_cmb_unlensed_solardipole_2Dorder2_LAT_MFL2_filtered_telescope_all_time_all_bmap.fits")
fname_cmb2 = os.path.join(rootdir, "chile_cmb_unlensed_solardipole_2Dorder5_LAT_MFL2_filtered_telescope_all_time_all_bmap.fits")
fname_cmb3 = os.path.join(rootdir, "chile_cmb_unlensed_solardipole_2Dorder10_LAT_MFL2_filtered_telescope_all_time_all_bmap.fits")
#fname_cmb4 = os.path.join(rootdir, "chile_cmb_unlensed_solardipole_2Dorder20_LAT_MFL2_filtered_telescope_all_time_all_bmap.fits")
fname_cmb4 = os.path.join(rootdir, "chile_cmb_unlensed_solardipole_2Dorder0_by_tube_LAT_MFL2_filtered_telescope_all_time_all_bmap.fits")
fname_cmb5 = os.path.join(rootdir, "chile_cmb_unlensed_solardipole_2Dorder1_by_tube_LAT_MFL2_filtered_telescope_all_time_all_bmap.fits")

fname_cmb_input = "/global/cscratch1/sd/zonca/cmbs4/map_based_simulations/202102_design_tool_input/4096/cmb_unlensed_solardipole/0000/cmbs4_cmb_unlensed_solardipole_uKCMB_LAT-MFL2_nside4096_0000.fits"

for fname_map in (
        fname_atmo00, fname_atmo0, fname_atmo1, fname_atmo2, fname_atmo3, fname_atmo4, fname_atmo5,
        fname_cmb00,fname_cmb0, fname_cmb1, fname_cmb2, fname_cmb3, fname_cmb4, fname_cmb5,
        fname_cmb_input):
    get_cl(fname_map, fname_hits)

"""
nrow, ncol = 2, 6
fig = plt.figure(figsize=[4 * ncol, 5 * nrow])
reso = 1
xsize = 800
cl_cmb_in = get_cl(fname_cmb_input)
for icol, (fname1, fname2, name) in enumerate([
        (fname_cmb00, fname_atmo00, "2Dorder=NO"),
        (fname_cmb0, fname_atmo0, "2Dorder=1"),
        (fname_cmb1, fname_atmo1, "2Dorder=2"),
        (fname_cmb2, fname_atmo2, "2Dorder=5"),
        (fname_cmb3, fname_atmo3, "2Dorder=10"),
        (fname_cmb4, fname_atmo3, "2Dorder=1 by tube"),
]):
    cl_cmb = get_cl(fname1)
    tf = cl_cmb_in[0] / cl_cmb[0]
    tf[lmax_tf:] = 1
    tf[:30] = 1
    m1 = hp.read_map(fname1)
    m2 = hp.read_map(fname2)
    for m in m1, m2,:
        alm = hp.map2alm(m, lmax=lmax, iter=0)
        alm = hp.almxfl(alm, tf ** .5)
        m[:] = hp.alm2map(alm, nside, lmax=lmax)

    hp.gnomview(m1, reso=reso, xsize=xsize, sub=[nrow, ncol, icol + 1], title="CMB " + name, min=-350, max=350)
    #hp.gnomview(m2, reso=reso, xsize=xsize, sub=[nrow, ncol, icol + 1 + ncol], title="Atmo " + name, min=-4e-3, max=4e-3)
    hp.gnomview(m2, reso=reso, xsize=xsize, sub=[nrow, ncol, icol + 1 + ncol], title="Atmo " + name, min=-0.005, max=0.005)
fig.savefig("stamps.deconv.png")
"""

nrow, ncol = 1, 2
fig = plt.figure(figsize=[6 * ncol, 6 * nrow])
ell = np.arange(lmax + 1)
norm = ell * (ell + 1) / (2 * np.pi) * 1e12

ax1 = fig.add_subplot(nrow, ncol, 1)
ax2 = fig.add_subplot(nrow, ncol, 2)

col = 0

cl_cmb_in = get_cl(fname_cmb_input)
ax2.loglog(ell[2:], (norm * cl_cmb_in[col])[2:], color="black", ls="-", label="CMB in")
    
cl_atmo = get_cl(fname_atmo00)
cl_cmb = get_cl(fname_cmb00)
tf = cl_cmb_in / cl_cmb
tf[:, lmax_tf:] = 1
#ax.loglog(ell[2:], (norm * cl_atmo[col])[2:], color="tab:blue", ls="--", lw=4)
ax1.loglog(ell[2:], (norm * cl_atmo[col])[2:] * tf[col][2:], color="tab:pink", ls="-", label="2Dorder=NO")
ax2.loglog(ell[2:], (norm * cl_cmb[col])[2:], color="tab:pink", ls="-", label="2Dorder=NO")
    
cl_atmo = get_cl(fname_atmo0)
cl_cmb = get_cl(fname_cmb0)
tf = cl_cmb_in / cl_cmb
tf[:, lmax_tf:] = 1
#ax.loglog(ell[2:], (norm * cl_atmo[col])[2:], color="tab:blue", ls="--", lw=4)
ax1.loglog(ell[2:], (norm * cl_atmo[col])[2:] * tf[col][2:], color="tab:blue", ls="-", label="2Dorder=1")
ax2.loglog(ell[2:], (norm * cl_cmb[col])[2:], color="tab:blue", ls="-", label="2Dorder=1")
    
#cl_atmo = get_cl(fname_atmo1)
#ax.loglog(ell[2:], (norm * cl_atmo[col])[2:], color="tab:orange", ls="--")
#ax.loglog(ell[2:], (norm * cl_atmo[col])[2:] * tf[col][2:], color="tab:orange", ls="-", label="2Xdet")

cl_atmo = get_cl(fname_atmo1)
cl_cmb = get_cl(fname_cmb1)
tf = cl_cmb_in / cl_cmb
tf[:, lmax_tf:] = 1
#ax.loglog(ell[2:], (norm * cl_atmo[col])[2:], color="tab:orange", ls="--")
ax1.loglog(ell[2:], (norm * cl_atmo[col])[2:] * tf[col][2:], color="tab:orange", ls="-", label="2Dorder=2")
ax2.loglog(ell[2:], (norm * cl_cmb[col])[2:], color="tab:orange", ls="-", label="2Dorder=2")

cl_atmo = get_cl(fname_atmo2)
cl_cmb = get_cl(fname_cmb2)
tf = cl_cmb_in / cl_cmb
tf[:, lmax_tf:] = 1
#ax.loglog(ell[2:], (norm * cl_atmo[col])[2:], color="tab:green", ls="--")
ax1.loglog(ell[2:], (norm * cl_atmo[col])[2:] * tf[col][2:], color="tab:green", ls="-", label="2Dorder=5")
ax2.loglog(ell[2:], (norm * cl_cmb[col])[2:], color="tab:green", ls="-", label="2Dorder=5")

cl_atmo = get_cl(fname_atmo3)
cl_cmb = get_cl(fname_cmb3)
tf = cl_cmb_in / cl_cmb
tf[:, lmax_tf:] = 1
#ax.loglog(ell[2:], (norm * cl_atmo[col])[2:], color="tab:green", ls="--")
ax1.loglog(ell[2:], (norm * cl_atmo[col])[2:] * tf[col][2:], color="tab:purple", ls="-", label="2Dorder=10")
ax2.loglog(ell[2:], (norm * cl_cmb[col])[2:], color="tab:purple", ls="-", label="2Dorder=10")

"""
cl_atmo = get_cl(fname_atmo4)
cl_cmb = get_cl(fname_cmb4)
tf = cl_cmb_in / cl_cmb
tf[:, lmax_tf:] = 1
#ax.loglog(ell[2:], (norm * cl_atmo[col])[2:], color="tab:green", ls="--")
ax1.loglog(ell[2:], (norm * cl_atmo[col])[2:] * tf[col][2:], color="tab:brown", ls="-", label="2Dorder=20")
ax2.loglog(ell[2:], (norm * cl_cmb[col])[2:], color="tab:brown", ls="-", label="2Dorder=20")
"""

cl_atmo = get_cl(fname_atmo4)
cl_cmb = get_cl(fname_cmb4)
tf = cl_cmb_in / cl_cmb
tf[:, lmax_tf:] = 1
#ax.loglog(ell[2:], (norm * cl_atmo[col])[2:], color="tab:green", ls="--")
ax1.loglog(ell[2:], (norm * cl_atmo[col])[2:] * tf[col][2:], color="tab:red", ls="-", label="2Dorder=0 by tube")
ax2.loglog(ell[2:], (norm * cl_cmb[col])[2:], color="tab:red", ls="-", label="2Dorder=0 by tube")

cl_atmo = get_cl(fname_atmo5)
cl_cmb = get_cl(fname_cmb5)
tf = cl_cmb_in / cl_cmb
tf[:, lmax_tf:] = 1
#ax.loglog(ell[2:], (norm * cl_atmo[col])[2:], color="tab:green", ls="--")
ax1.loglog(ell[2:], (norm * cl_atmo[col])[2:] * tf[col][2:], color="tab:brown", ls="-", label="2Dorder=1 by tube")
ax2.loglog(ell[2:], (norm * cl_cmb[col])[2:], color="tab:brown", ls="-", label="2Dorder=1 by tube")


for ax in ax1, ax2:
    ax.set_xlabel("Multipole, $\ell$")
    ax.set_ylabel("D$_\ell$")
    comp = ["TT", "EE", "BB"][col]
ax1.set_title(f"Atmosphere {comp}")
ax2.set_title(f"CMB {comp}")

ax.legend(loc="best")
fig.savefig("chlat_atmo_study.png")
