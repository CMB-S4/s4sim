import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np

import toast.qarray as qa

RADIUS = 10  # 15.0  # degrees
THROW = 40  # degrees
SCANTIME = 60  # minutes
XAXIS, YAXIS, ZAXIS = np.eye(3)

DEC_MIN = -70
DEC_MAX = 30

targets = {
    "South_Deep" : (30, -45, 1, 90),  # Ra, Dec, base weight, max_radius
    "North_Deep" : (165, 5, 1e10, 60),  # Ra, Dec, base weight, max_radius
}

fname_fg = "/home/reijo/work/npipe6/npipe6v20_353_map.fits"
fname_fg_smooth = "smooth_npipe6v20_353_map.fits"
if not os.path.isfile(fname_fg_smooth):
    print(f"Loading {fname_fg}")
    fg = hp.read_map(fname_fg, None)
    lmax = 1024
    alm = hp.map2alm(fg, lmax=lmax, iter=0)
    rot = hp.Rotator(coord="GC")
    alm = rot.rotate_alm(alm)
    sfg = hp.alm2map(alm, 512, fwhm=np.radians(3))
    print(f"Writing {fname_fg_smooth}")
    hp.write_map(fname_fg_smooth, sfg)
else:
    print(f"Loading {fname_fg_smooth}")
    sfg = hp.read_map(fname_fg_smooth, None) * 1e6

p = np.sqrt(sfg[1]**2 + sfg[2]**2)
plimit = 100
hp.mollview(p, min=0, max=plimit, cmap="inferno")
nside = hp.get_nside(p)

fname_targets = "max-depth-targets.alt01.txt"
f = open(fname_targets, "w")

for primary in targets:
    ra, dec, base_weight, max_radius = targets[primary]

    hp.projplot([ra], [dec], 'o', color="red", lonlat=True)
    f.write(f"--patch\n{primary},MAX-DEPTH,{base_weight},{ra:.1f},{dec:.1f},{RADIUS},{THROW},{SCANTIME}\n")

    quat = qa.from_lonlat_angles(np.radians(ra), np.radians(dec), 0)
    for radius in range(RADIUS, max_radius + 1, RADIUS):
        alpha = 1 - radius / 100
        r = np.radians(radius)
        vec_in = np.array([np.sin(r), 0, np.cos(r)])
        lons, lats = [], []
        n = int(2 * np.pi * radius / RADIUS)
        for angle in np.linspace(0, 360, n, endpoint=False):
            zrot = qa.from_axisangle(ZAXIS, np.radians(angle))
            vec_out = qa.rotate(quat, qa.rotate(zrot, vec_in))
            pix = hp.vec2pix(nside, *vec_out)
            if p[pix] > plimit:
                continue
            lon, lat = hp.vec2dir(vec_out, lonlat=True)
            if lat < DEC_MIN or lat > DEC_MAX:
                continue
            lons.append(lon)
            lats.append(lat)

        hp.projplot(lons, lats, 'o', color="orange", alpha=alpha, lonlat=True)

        weight = base_weight * 10**radius
        for i, (lon, lat) in enumerate(zip(lons, lats)):
            name = f"{primary}-{radius}-{i:03}"
            f.write(f"--patch\n{name},MAX-DEPTH,{weight:.1e},{lon:.1f},{lat:.1f},{RADIUS},{THROW},{SCANTIME}\n")
fname_plot = "max-depth-targets.png"
plt.savefig(fname_plot)
print(f"Wrote {fname_plot}")

f.close()
print(f"Wrote {fname_targets}")
