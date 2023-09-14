import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse

from toast.pixels_io_healpix import read_healpix


rootdir = "/global/cfs/cdirs/cmbs4/dc/dc0/staging"
nside = 512
npix = 12 * nside**2
lmax = 2 * nside

nsplit = 32
if len(sys.argv) > 1:
    band = int(sys.argv[1])
else:
    band = 85
if len(sys.argv) > 2:
    band_obsmat = int(sys.argv[2])
else:
    band_obsmat = band
fname_input = f"../cmb_sim/input_maps/cmb.spsat.f{band:03}.h5"
fname_matrix = f"{rootdir}/obsmat/obsmats/full/f{band_obsmat:03}/obsmat_f{band_obsmat:03}.npz"
fname_cov = f"{rootdir}/noise_sim/outputs_rk/coadd/spsat/coadd_spsat_f{band:03}_001of001_cov.fits"

print(f"Loading {fname_matrix}")
obs_matrix = scipy.sparse.load_npz(fname_matrix)

print(f"Loading {fname_input}")
input_map = read_healpix(fname_input, None, nest=True)

print("Applying matrix")
test_map = obs_matrix.dot(input_map.ravel()).reshape([3, -1])
test_map = hp.reorder(test_map, n2r=True)

print(f"Loading {fname_cov}")
w = hp.read_map(fname_cov)
good = w > 0
sorted_w = np.sort(w[good])
lim = sorted_w[int(npix * 0.03)]
mask = np.ones(npix, dtype=bool)
mask[w == 0] = False
mask[w > lim] = False

input_map = hp.reorder(input_map, n2r=True)
input_map[0] = hp.remove_dipole(input_map[0] * mask, bad=0)
cl_input = hp.anafast(input_map * mask, lmax=lmax, iter=0)
cl_test = hp.anafast(test_map * mask, lmax=lmax, iter=0)

cl_output = []
cl_diff = []
for isplit in range(1, nsplit + 2):
    if isplit > nsplit:
        # Total mission result
        fname_output = f"{rootdir}/multimap_sim/outputs_rk/coadd/spsat/coadd_spsat_f{band:03}_unlensed_cmb_001of001_map.fits"
    else:
        fname_output = f"{rootdir}/multimap_sim/outputs_rk/coadd/spsat/coadd_spsat_f{band:03}_unlensed_cmb_{isplit:03}of{nsplit:03}_map.fits"
    print(f"Loading {fname_output}")
    output_map = read_healpix(fname_output, None)

    cl_output.append(hp.anafast(output_map * mask, lmax=lmax, iter=0))
    cl_diff.append(hp.anafast((test_map - output_map) * mask, lmax=lmax, iter=0))

ell = np.arange(lmax + 1)
nrow, ncol = 1, 3
fig = plt.figure(figsize=[4 * ncol, 4 * nrow])
for i in range(3):
    ax = fig.add_subplot(nrow, ncol, 1 + i)
    ax.loglog(ell[2:], cl_input[i][2:] * 1e12, color="black", label="Input")
    ax.loglog(ell[2:], cl_test[i][2:] * 1e12, color="tab:blue", label="Obsmat X Input")
    label1 = f"DC0, nsplit={nsplit}"
    label2 = f"Obsmat - DC0, nsplit={nsplit}"
    for isplit in range(nsplit):
        ax.loglog(ell[2:], cl_output[isplit][i][2:] * 1e12, color="tab:orange", label=label1)
        ax.loglog(ell[2:], cl_diff[isplit][i][2:] * 1e12, color="tab:green", label=label2)
        label1 = None
        label2 = None
    ax.loglog(ell[2:], cl_output[-1][i][2:] * 1e12, color="tab:red", label="DC0, full")
    ax.loglog(ell[2:], cl_diff[-1][i][2:] * 1e12, color="tab:purple", label="Obsmat - DC0, full")
    comp = ["TT", "EE", "BB"][i]
    ax.set_title(f"{comp} {band}GHz")
    ax.set_xlabel("Multipole, $\ell$")
    ax.set_ylabel("C$\ell$ [$\mu$K$^2$]")
ax.legend(loc="best")

fig.tight_layout()
if band == band_obsmat:
    fname_plot = f"obsmat_cl_comparison_f{band:03}_nsplit{nsplit:03}.png"
else:
    fname_plot = f"obsmat_cl_comparison_f{band:03}_f{band_obsmat:03}_nsplit{nsplit:03}.png"
fig.savefig(fname_plot)
print(f"Wrote {fname_plot}")
