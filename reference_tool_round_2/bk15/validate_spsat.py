import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np

import pylab

import requirements as req

fiducial_ellnorm = req.fiducial_ellnorm

telescope = "SAT"
bands = "LFS1", "LFS2", "MFLS1", "MFLS2", "MFHS1", "MFHS2", "HFS1", "HFS2"
nband = len(bands)
site = "pole"
nside = 512
lmax = 2 * nside
flavor = "noise_atmo_7splits" # cmb_r0  cmb_tensor_only_r3e-3  foregrounds  noise_atmo_7splits

split = 1
nsplits = 1

bk15 = np.genfromtxt("BK15_Nell", skip_header=13).T
bk15_ellnorm = bk15[0] * (bk15[0] + 1) / (2 * np.pi)
bl_bk15_95 = np.exp(-bk15[0] * (bk15[0] + 1) * np.radians(0.30) ** 2)
bl_bk15_150 = np.exp(-bk15[0] * (bk15[0] + 1) * np.radians(0.21) ** 2)
bl_bk15_220 = np.exp(bk15[0] * (bk15[0] + 1) * np.radians(0.14) ** 2)
bk15[1:4] *= bl_bk15_95 / bk15_ellnorm
bk15[4:7] *= bl_bk15_150 / bk15_ellnorm
bk15[7:10] *= bl_bk15_220 / bk15_ellnorm


def map2cl(m, hits):
    good = hits > 0
    ngood = np.sum(good)
    sorted_hits = np.sort(hits[good])
    hit_lim = sorted_hits[np.int(ngood * .5)]
    mask = hits > hit_lim
    fsky = np.sum(mask) / mask.size
    cl = hp.anafast(m * mask, lmax=lmax, iter=0) / fsky

    return cl


def get_cl(fname_cl, fname_map, fname_hits):
    if os.path.isfile(fname_cl):
        cl = hp.read_cl(fname_cl)
    else:
        hits = hp.read_map(path_hits)
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
        hits = hp.read_map(fname_hits)
        cl_in = map2cl(inmap, hits)
        cl_out = map2cl(outmap, hits)
        tf = cl_out / cl_in
        hp.write_cl(fname_tf, tf)
    return tf


rootdir = "/global/cscratch1/sd/zonca/cmbs4/map_based_simulations/202102_design_tool_run"
rootdir_bk15 = "/global/cscratch1/sd/zonca/cmbs4/map_based_simulations/202103_bicepkeck/noise_atmo"

nrow, ncol = nband, 3
fig1 = plt.figure(figsize=[6 * ncol, 6 * 1])
ax1 = fig1.add_subplot(1, 3, 1)
ax2 = fig1.add_subplot(1, 3, 2)
ax3 = fig1.add_subplot(1, 3, 3)
fig2 = plt.figure(figsize=[6 * ncol, 6 * nrow])
ell = np.arange(lmax + 1)
# ellnorm = ell * (ell + 1) / (2 * np.pi) * 1e12
ellnorm = 1e12
iplot = 0

for band in bands:
    path_hits = os.path.join(
        rootdir,
        "noise_atmo_7splits",
        f"{telescope}-{band}_{site}/cmbs4_hitmap_{telescope}-{band}_{site}_nside{nside}_{split}_of_{nsplits}.fits",
    )
    # Transfer function
    path_tf = f"tf_{telescope}_{band}_{site}.fits"
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
    path_cl = f"cl_{telescope}_{band}_{site}_{flavor}.fits"
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
    fwhm = req.Pole_SAT[freq][0]
    bl = req.get_bl(fwhm, ell)
    nltt = req.NlTT_Pole_SAT[freq]
    nlee = req.NlEE_Pole_SAT[freq]

    iplot += 1
    ax = fig2.add_subplot(nrow, ncol, iplot)
    ax.set_title(f"TT {band} / {freq}GHz")
    ax.set_xlabel("Multipole, $\ell$")
    ax.set_ylabel("D$\ell$ [$\mu$K$^2$]")
    ax.loglog(req.fiducial_ell, req.fiducial_TT / fiducial_ellnorm, "k", label="CMB")
    ax.loglog(req.ells, nltt / req.ellnorm, label="requirement")
    ax.loglog(ell, ellnorm * cl[0] * bl, label=f"Sim")
    if band == "MFLS1" or band == "MFHS1":
        ax.loglog(bk15[0], bk15[1], label="BK15 95GHz")
    elif band == "MFLS2" or band == "MFHS2":
        ax.loglog(bk15[0], bk15[4], label="BK15 150GHz")
    elif band == "HFS1":
        ax.loglog(bk15[0], bk15[7], label="BK15 220GHz")
    #ax.set_xscale("linear")
    #ax.set_yscale("log")
    ax.set_xlim([10, 1000])
    ax.set_ylim([1e-8, 1e6])

    iplot += 1
    ax = fig2.add_subplot(nrow, ncol, iplot)
    ax.set_title(f"EE {band} / {freq}GHz")
    ax.set_xlabel("Multipole, $\ell$")
    ax.set_ylabel("D$\ell$ [$\mu$K$^2$]")
    ax.loglog(req.fiducial_ell, req.fiducial_EE / fiducial_ellnorm, "k", label="CMB")
    ax.loglog(req.ells, nlee / req.ellnorm, label="requirement")
    ax.loglog(ell, ellnorm * cl[1] * bl, label=f"Sim")
    if band == "MFLS1" or band == "MFHS1":
        ax.loglog(bk15[0], bk15[2], label="BK15 95GHz")
    elif band == "MFLS2" or band == "MFHS2":
        ax.loglog(bk15[0], bk15[5], label="BK15 150GHz")
    elif band == "HFS1":
        ax.loglog(bk15[0], bk15[8], label="BK15 220GHz")
    #ax.set_xscale("linear")
    #ax.set_yscale("log")
    ax.set_xlim([10, 1000])
    ax.set_ylim([1e-8, 1e0])
    
    iplot += 1
    ax = fig2.add_subplot(nrow, ncol, iplot)
    ax.set_title(f"BB {band} / {freq}GHz")
    ax.set_xlabel("Multipole, $\ell$")
    ax.set_ylabel("D$\ell$ [$\mu$K$^2$]")
    ax.loglog(req.fiducial_ell, req.fiducial_BB / fiducial_ellnorm, "k", label="CMB")
    ax.loglog(req.ells, nlee / req.ellnorm, label="requirement")
    ax.loglog(ell, ellnorm * cl[2] * bl, label=f"Sim")
    if band == "MFLS1" or band == "MFHS1":
        ax.loglog(bk15[0], bk15[3], label="BK15 95GHz")
    elif band == "MFLS2" or band == "MFHS2":
        ax.loglog(bk15[0], bk15[6], label="BK15 150GHz")
    elif band == "HFS1":
        ax.loglog(bk15[0], bk15[9], label="BK15 220GHz")
    ax.legend(loc="best")
    #ax.set_xscale("linear")
    #ax.set_yscale("log")
    ax.set_xlim([10, 1000])
    ax.set_ylim([1e-8, 1e0])

for ax in [ax1, ax2, ax3]:
    ax.set_xlabel(r"Multipole, $\ell$")
    ax.set_ylabel("Transfer function")
    ax.set_xlim([1, lmax + 1])
    ax.set_ylim([-0.1, 1.1])
    ax.axhline(1.0, color="k", linestyle="--")
    ax.set_xscale("log")
ax3.legend(loc="best")
fig1.savefig("pole_sat_tf.png")

fig2.savefig(f"pole_sat_validation.{flavor}.png")
plt.show()
