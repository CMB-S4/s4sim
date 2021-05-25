import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
import pylab

import requirements as req

#pylab.rc('text', usetex=True)
#pylab.rc('font', family='serif')

telescope = "LAT"
#bands = "LFL1", "LFL2", "MFL1", "MFL2", "HFL1", "HFL2"
bands = "MFL1", "MFL2"
nband = len(bands)
site = "chile"
nside = 4096
lmax = 2 * nside
flavor = "noise_atmo_7splits" # cmb_r0  cmb_tensor_only_r3e-3  foregrounds  noise_atmo_7splits
#flavor = "noise" # cmb_r0  cmb_tensor_only_r3e-3  foregrounds  noise_atmo_7splits

split = 1
nsplits = 1

# Mock the ACT processing mask.  This is 525 sq.deg (0.0127 fsky)

npix = 12 * nside ** 2
pix = np.arange(npix)
vec = hp.pix2vec(nside, pix)
lon, lat = hp.vec2dir(vec, lonlat=True)
mask = np.logical_and(lon > -4, lon < 40) * np.logical_and(lat > -8, lat < 4)

# Deep56 is 565 sq.deg (0.0137 fsky) of which 340 sq.deg (0.00824 fsky) is usable for power spectrum estimation
# ell pa1(150GHz) pa2(150GHz) pa3(150GHz) pa3(98GHz) pa3(98x150GHz)
act_tt = np.genfromtxt("deep56_TT_Nl_out_210317.txt", skip_header=1).T
act_ee = np.genfromtxt("deep56_EE_Nl_out_210317.txt", skip_header=1).T


def map2cl(m, hits):
    """
    good = hits > 0
    ngood = np.sum(good)
    sorted_hits = np.sort(hits[good])
    
    hit_lim = sorted_hits[np.int(ngood * .01)]
    mask = hits > hit_lim
    """
    fsky = np.sum(mask) / mask.size
    m[0] -= np.mean(m[0][mask])
    cl = hp.anafast(m * mask, lmax=lmax, iter=0) / fsky
    
    return cl


def get_cl(fname_cl, fname_map, fname_hits):
    if os.path.isfile(fname_cl):
        cl = hp.read_cl(fname_cl)
    else:
        hits = None # hp.read_map(path_hits)
        m = hp.read_map(fname_map, None)
        cl = map2cl(m, hits)
        hp.write_cl(fname_cl, cl)
    return cl


def get_tf(fname_tf, fname_cmb_unlensed, fname_cmb_lensing, fname_output, fname_hits):
    if os.path.isfile(fname_tf):
        tf = hp.read_cl(fname_tf)
    else:
        inmap = hp.read_map(fname_cmb_unlensed, None) + hp.read_map(fname_cmb_lensing, None)
        inmap *= 1e-6  # into K_CMB
        inmap[0] = hp.remove_dipole(inmap[0])
        outmap = hp.read_map(fname_output, None)
        hits = None # hp.read_map(fname_hits)
        cl_in = map2cl(inmap, hits)
        cl_out = map2cl(outmap, hits)
        tf = cl_out / cl_in
        hp.write_cl(fname_tf, tf)
    return tf


rootdir = "/global/cscratch1/sd/zonca/cmbs4/map_based_simulations/202102_design_tool_run"

nrow, ncol = nband, 3
fig1 = plt.figure(figsize=[6 * ncol, 6 * 1])
ax1 = fig1.add_subplot(1, 3, 1)
ax2 = fig1.add_subplot(1, 3, 2)
ax3 = fig1.add_subplot(1, 3, 3)
fig2 = plt.figure(figsize=[6 * ncol, 6 * nrow])
ell = np.arange(lmax + 1)
ellnorm = ell * (ell + 1) / (2 * np.pi) * 1e12
iplot = 0

for band in bands:
    path_hits = os.path.join(
        rootdir,
        "noise_atmo_7splits",
        f"{telescope}-{band}_{site}/cmbs4_hitmap_{telescope}-{band}_{site}_nside{nside}_{split}_of_{nsplits}.fits",
    )
    # Transfer function
    path_tf = f"tf_{telescope}_{band}_{site}.deep56.fits"
    path_cmb_unlensed = os.path.join(
        "/global/cscratch1/sd/zonca/cmbs4/map_based_simulations/202102_design_tool_input",
        f"{nside}/cmb_unlensed_solardipole/0000/"
        f"cmbs4_cmb_unlensed_solardipole_uKCMB_{telescope}-{band}_nside{nside}_0000.fits",
    )
    path_cmb_lensing = os.path.join(
        "/global/cscratch1/sd/zonca/cmbs4/map_based_simulations/202102_design_tool_input",
        f"{nside}/cmb_lensing_signal/0000/"
        f"cmbs4_cmb_lensing_signal_uKCMB_{telescope}-{band}_nside{nside}_0000.fits",
    )
    path_cmb_output = os.path.join(
        rootdir,
        "cmb_r0",
        f"{telescope}-{band}_{site}/cmbs4_KCMB_{telescope}-{band}_{site}_nside{nside}_1_of_1.fits",
    )
    tf = get_tf(path_tf, path_cmb_unlensed, path_cmb_lensing, path_cmb_output, path_hits)
    # N_ell
    path_cl = f"cl_{telescope}_{band}_{site}_{flavor}.deep56.fits"
    path_noise_map = os.path.join(
        rootdir,
        flavor,
        f"{telescope}-{band}_{site}/"
        f"cmbs4_KCMB_{telescope}-{band}_{site}_nside{nside}_{split}_of_{nsplits}.fits",
    )
    cl = get_cl(path_cl, path_noise_map, path_hits) / tf

    ax1.plot(ell[2:], tf[0][2:], label=band)
    ax2.plot(ell[2:], tf[1][2:], label=band)
    ax3.plot(ell[2:], tf[2][2:], label=band)

    freq = req.band2freq[band]
    fwhm = req.Chile_LAT[freq][0]
    bl = req.get_bl(fwhm, ell)
    nltt = req.NlTT_Chile_LAT[freq]
    nlee = req.NlEE_Chile_LAT[freq]

    act_obs_eff = 0.1  # fraction of dedicated to deep56 in 2015, includes cutting daytime observations
    scale_t = 7 / 1 * (1 / 68) / act_obs_eff
    scale_p = scale_t * 100  # last factor accounts for differences in FP size

    iplot += 1
    ax = fig2.add_subplot(nrow, ncol, iplot)
    ax.set_title(f"TT {band} / {freq}GHz")
    ax.set_xlabel("Multipole, $\ell$")
    ax.set_ylabel("D$\ell$ [$\mu$K$^2$]")
    ax.loglog(req.fiducial_ell, req.fiducial_TT, "k", label="CMB")
    #ax.loglog(req.ells, nltt, label="requirement")
    ax.loglog(ell, ellnorm * cl[0] * bl * scale_t, label=f"Sim")
    if band == "MFL1":
        ax.loglog(act_tt[0], act_tt[4], label="ACT PA3")
    elif band == "MFL2":
        ax.loglog(act_tt[0], act_tt[1], label="ACT PA1")
        ax.loglog(act_tt[0], act_tt[2], label="ACT PA2")
        ax.loglog(act_tt[0], act_tt[3], label="ACT PA3")
    #ax.set_xscale("linear")
    #ax.set_yscale("log")
    ax.set_xlim([20, 8000])
    ax.set_ylim([1e-1, 1e7])

    iplot += 1
    ax = fig2.add_subplot(nrow, ncol, iplot)
    ax.set_title(f"EE {band} / {freq}GHz")
    ax.set_xlabel("Multipole, $\ell$")
    ax.set_ylabel("C$\ell$ [$\mu$K$^2$]")
    ax.loglog(req.fiducial_ell, req.fiducial_EE, "k", label="CMB")
    #ax.loglog(req.ells, nlee, label="requirement")
    ax.loglog(ell, ellnorm * cl[1] * bl * scale_p, label=f"Sim")
    if band == "MFL1":
        ax.loglog(act_ee[0], act_ee[4], label="ACT PA3")
    elif band == "MFL2":
        ax.loglog(act_ee[0], act_ee[1], label="ACT PA1")
        ax.loglog(act_ee[0], act_ee[2], label="ACT PA2")
        ax.loglog(act_ee[0], act_ee[3], label="ACT PA3")
    #ax.set_xscale("linear")
    #ax.set_yscale("log")
    ax.set_xlim([20, 8000])
    ax.set_ylim([1e-4, 1e5])

    iplot += 1
    ax = fig2.add_subplot(nrow, ncol, iplot)
    ax.set_title(f"BB {band} / {freq}GHz")
    ax.set_xlabel("Multipole, $\ell$")
    ax.set_ylabel("C$\ell$ [$\mu$K$^2$]")
    ax.loglog(req.fiducial_ell, req.fiducial_BB, "k", label="CMB")
    #ax.loglog(req.ells, nlee, label="requirement")
    ax.loglog(ell, ellnorm * cl[2] * bl * scale_p, label=f"Sim")
    if band == "MFL1":
        ax.loglog(act_ee[0], act_ee[4], label="ACT PA3")
    elif band == "MFL2":
        ax.loglog(act_ee[0], act_ee[1], label="ACT PA1")
        ax.loglog(act_ee[0], act_ee[2], label="ACT PA2")
        ax.loglog(act_ee[0], act_ee[3], label="ACT PA3")
    ax.legend(loc="best")
    #ax.set_xscale("linear")
    #ax.set_yscale("log")
    ax.set_xlim([20, 8000])
    ax.set_ylim([1e-4, 1e5])

for ax in [ax1, ax2, ax3]:
    ax.set_xlabel(r"Multipole, $\ell$")
    ax.set_ylabel("Transfer function")
    ax.set_xlim([1, lmax + 1])
    ax.set_ylim([-0.1, 1.1])
    ax.axhline(1.0, color="k", linestyle="--")
    ax.set_xscale("log")
ax3.legend(loc="best")
fig1.savefig("chile_lat_tf.deep56.png")

fig2.savefig(f"chile_lat_validation.{flavor}.deep56.png")
plt.show()
