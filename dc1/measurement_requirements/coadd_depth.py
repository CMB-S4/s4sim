import os
import sys

import healpy as hp
import numpy as np


fsky = 0.60
all_mr = {
    30 : (21.8, 30.8),
    40 : (12.4, 17.6),
    90 : (2.0, 2.9),
    150 : (2.0, 2.8),
    220 : (6.9, 9.8),
    280 : (16.7, 23.6),
}

rootdir = "/global/cfs/cdirs/cmbs4/dc/dc1/staging/noise_sim/outputs_rk/coadd/LAT0_CHLAT"
for freq in 30, 40, 90, 150, 220, 280:
    band = f"f{freq:03}"
    print(f"\n{band}")
    fname_cov = os.path.join(rootdir, f"coadd_LAT0_CHLAT_{band}_001of001_cov.fits")
    print(f"Loading {fname_cov}")
    ii, qq, uu = hp.read_map(fname_cov, [0, 3, 5])
    nside = hp.get_nside(ii)
    pixarea = hp.nside2pixarea(nside, degrees=True)
    tdepth = np.sqrt(ii * pixarea) * 1e6 * 60
    qdepth = np.sqrt(qq * pixarea) * 1e6 * 60
    udepth = np.sqrt(uu * pixarea) * 1e6 * 60
    pdepth = np.fmax(qdepth, udepth)
    mr_t, mr_p = all_mr[freq]
    limits = []
    for mr, depth in zip([mr_t, mr_p], [tdepth, pdepth]):
        depth[depth == 0] = np.inf
        depth = np.sort(depth)
        npix = depth.size
        i = int(fsky * npix)
        limits.append(depth[i])
        print(f"Depth at fsky = {fsky} is {depth[i]:.3f} uK.arcmin.  ", end="")
        i = np.searchsorted(depth, mr)
        print(f"fsky that meets depth < {mr} uK.arcmin is {i / npix:.03f}")
    fname_map = os.path.join(rootdir, f"coadd_LAT0_CHLAT_{band}_001of001_map.fits")
    print(f"Loading {fname_map}")
    m = hp.read_map(fname_map, None)
    best = tdepth < limits[0]
    rms_i = np.std(m[0][best])
    best = pdepth < limits[1]
    rms_q = np.std(m[1][best])
    rms_u = np.std(m[2][best])
    scale = np.sqrt(pixarea) * 1e6 * 60
    print(f"Simulated depth (best {fsky*100}%):")
    for stokes, rms in zip("IQU", [rms_i, rms_q, rms_u]):
        print(f"{stokes} = {rms * scale:.03f} uK.arcmin")
