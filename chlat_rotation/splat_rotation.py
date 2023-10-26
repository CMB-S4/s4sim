import os
import pickle
import sys

import astropy.units as u
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np

import toast
import toast.io
import toast.qarray as qa


tele = "SPLAT"

XAXIS, YAXIS, ZAXIS = np.eye(3)

if len(sys.argv) > 1:
    freq = int(sys.argv[1])
else:
    #freq = 30
    #freq = 90
    freq = 220
band = f"f{freq:03}"
fname = f"focalplane_LAT2_SPLAT_{band}.h5"

comm, procs, rank = toast.get_world()

fp = toast.instrument.Focalplane()
print(f"Loading {fname}")
with toast.io.H5File(fname, "r", comm=comm) as f:
    fp.load_hdf5(f.handle, comm=comm)

quats =  fp.detector_data["quat"]

def count_empty(hist):
    ind = np.ones(hist.size, dtype=bool)
    istart = 0
    while istart < hist.size:
        if hist[istart] != 0:
            break
        istart += 1
    istop = hist.size
    while istop > 0:
        if hist[istop - 1] != 0:
            break
        istop -= 1
    frac = np.sum(hist[istart:istop] != 0) / (istop - istart)
    rms = np.std(hist[istart:istop])
    return frac, rms

fp_radius = None
R = None
fwhm = np.median(fp.detector_data["fwhm"]).to_value(u.deg)
if freq < 50:
    wbin = (2.0 * u.arcmin).to_value(u.deg)
elif freq < 200:
    wbin = (1.0 * u.arcmin).to_value(u.deg)
else:
    wbin = (0.5 * u.arcmin).to_value(u.deg)
nbin = None
target = None
ymax = None

#if freq > 50 and freq < 200:
if freq < 200:
    # The polarization angles are 15 degrees apart starting at 7.5 deg.
    # 7.5 + [0, 15, 30, 45, 60, 75]
    npol = 6
    psi_pol = fp.detector_data["pol_ang"]
    pol_type = (np.around((((psi_pol - 7.5) % 90) / 15)) % npol).astype(int)
    pol_ang = 7.5 + np.arange(npol) * 15
elif freq > 200:
    npol = 2
    psi_pol = fp.detector_data["pol_ang"]
    pol_type = (np.around(((psi_pol % 90) / 45)) % npol).astype(int)
    pol_ang = np.arange(npol) * 45
#else:
#    raise RuntimeError(f"Unknown freq = {freq}")

#for angle in 0, 1, 2, 3, 4, 5, 10, 15, 20, 22.5, 30, 45, 60, 67.5, 90:
#angles = np.arange(91)
#angles = np.linspace(0, 90, 901)
angles = np.linspace(0, 90, 9001)
nangle = angles.size
frac_vec_lat = np.zeros(nangle)
rms_vec_lat = np.zeros(nangle)
frac_vec_lon = np.zeros(nangle)
rms_vec_lon = np.zeros(nangle)
frac_vec_pol = np.zeros([npol, nangle])
rms_vec_pol = np.zeros([npol, nangle])

for iangle, angle in enumerate(angles):
    if np.abs(np.round(angle) - angle) < 1e-3:
        angle = np.round(angle)
        do_plot = True
    else:
        do_plot = False
    print(f"Converting quats to angles, angle = {angle} deg")
    yrot = qa.rotation(YAXIS, np.pi / 2)
    zrot = qa.rotation(ZAXIS, np.radians(angle))
    quats2 = qa.mult(yrot, qa.mult(zrot, quats))
    theta, phi, psi = qa.to_iso_angles(quats2)
    lon = np.degrees(phi)
    lat = 90 - np.degrees(theta)

    if fp_radius is None:
        # fp_radius = max(np.amax(np.abs(lon)), np.amax(np.abs(lat)))
        fp_radius = np.amax(np.sqrt(lon**2 + lat**2))
        R = int(fp_radius) + 1
        nbin = int(2 * R / wbin)
        target = len(lon) / (2 * fp_radius / wbin)

    lon_hist, lon_edges = np.histogram(lon, bins=nbin, range=[-R, R])
    lat_hist, lat_edges = np.histogram(lat, bins=nbin, range=[-R, R])
    lat_hists = []
    for pol in range(npol):
        ind = pol_type == pol
        lat_hists.append(np.histogram(lat[ind], bins=nbin, range=[-R, R])[0])
    # bin_centers = 0.5 * (lon_edges[:-1] + lon_edges[1:])
    if ymax is None:
        ymax = (np.amax(lat_hist) // 10 + 1) * 10

    if do_plot:
        nrow, ncol = 1, 3
        fig = plt.figure(figsize=[8 * ncol, 8 * nrow])
        fig.suptitle(f"Rotation = {angle} deg, wbin = {wbin * 60:.2f} arcmin")
        ax = fig.add_subplot(nrow, ncol, 1)
        ax.plot(lon, lat, '.')
        ax.set_xlabel("In-scan [deg]")
        ax.set_ylabel("Cross-scan [deg]")
        fov = 10.0
        ax.set_xlim([-fov/2, fov/2])
        ax.set_ylim([-fov/2, fov/2])
        ax.grid()
        ax.set_aspect("equal")

    frac, rms = count_empty(lon_hist)
    frac_vec_lon[iangle] = frac
    rms_vec_lon[iangle] = rms
    if do_plot:
        ax = fig.add_subplot(nrow, ncol, 2)
        ax.set_title(f"In-scan {frac:.3f}, $\sigma$={rms:6.1f}")
        ax.stairs(lon_hist, lon_edges)
        ax.axhline(target, color="k", linestyle="--")
        ax.set_ylim([0, ymax])
        ax.set_xlabel("In-scan [deg]")
        ax.set_ylabel("Ndetectors")

    frac, rms = count_empty(lat_hist)
    frac_vec_lat[iangle] = frac
    rms_vec_lat[iangle] = rms
    if do_plot:
        ax = fig.add_subplot(nrow, ncol, 3)
        ax.set_title(f"Cross-scan {frac:.3f}, $\sigma$={rms:6.1f}")
        ax.stairs(lat_hist, lat_edges)

    for ptype, hist in enumerate(lat_hists):
        frac, rms = count_empty(hist)
        frac_vec_pol[ptype, iangle] = frac
        rms_vec_pol[ptype, iangle] = rms
        if do_plot:
            ax.stairs(hist, lat_edges, label=f"{frac:.3f}, $\sigma$={rms:5.1f}")

    if do_plot:
        ax.axhline(target, color="k", linestyle="--")
        ax.set_ylim([0, ymax])
        ax.legend(loc="upper right")
        ax.set_xlabel("Cross-scan [deg]")
        ax.set_ylabel("Ndetectors")

        fname_plot = f"histograms_{tele}_{band}_{angle:0>5.1f}deg.png"
        fig.tight_layout()
        fig.savefig(fname_plot)
        print(f"Saved plot to {fname_plot}")
        plt.close()

fname_data = f"data_{tele}_{band}.pck"
with open(fname_data, "wb") as f:
    pickle.dump([tele, band, pol_ang, angles, frac_vec_lon, frac_vec_lat, frac_vec_pol, rms_vec_lon, rms_vec_lat, rms_vec_pol], f)

nrow, ncol = 2, 2
fig = plt.figure(figsize=[8 * ncol, 8 * nrow])
#fig.suptitle(f"Rotation = {angle} deg")

ax1 = fig.add_subplot(nrow, ncol, 1)
ax1.plot(angles, frac_vec_lon, '-', label="Co-scan fraction")
ax1.plot(angles, frac_vec_lat, '-', label="Cross-scan fraction")
ax1.legend(loc="best")

ax2 = fig.add_subplot(nrow, ncol, 2)
ax2.plot(angles, rms_vec_lon, '-', label="Co-scan RMS")
ax2.plot(angles, rms_vec_lat, '-', label="Cross-scan RMS")
ax2.legend(loc="best")

ax3 = fig.add_subplot(nrow, ncol, 3)
ax3.set_title("Cross-scan fraction by pol type")
for ipol in range(npol):
    ang = pol_ang[ipol]
    ax3.plot(angles, frac_vec_pol[ipol], '-', label=f"{ang} deg")
worst_frac = np.amin(frac_vec_pol, 0)
ax3.plot(angles, worst_frac, 'k--', label="worst")
#ax3.legend(loc="best")

ax4 = fig.add_subplot(nrow, ncol, 4)
ax4.set_title("Cross-scan RMS by pol type")
for ipol in range(npol):
    ang = pol_ang[ipol]
    ax4.plot(angles, rms_vec_pol[ipol], '-', label=f"{ang} deg")
worst_rms = np.amax(rms_vec_pol, 0)
ax4.plot(angles, worst_rms, 'k--', label="Worst")
ax4.legend(loc="best")

ibest = np.argmin(worst_rms)
best_angle = angles[ibest]
for ax in [ax1, ax2, ax3, ax4]:
    ax.axvline(best_angle, color="r", linestyle="--")

fig.suptitle(
    f"band = {band}, Best angle = {best_angle:.3f} deg, Bin width = {wbin * 60:.2f} arcmin"
)

fig.tight_layout()
fname = f"fraction_and_rms_{tele}_{band}.png"
fig.savefig(fname)
print(f"Saved plot in {fname}")
