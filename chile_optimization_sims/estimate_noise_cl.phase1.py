# Study the performance-based scaling applied in Phase-1 simulations

# Measure the per-pixel noise variance from the MC noise maps

import glob
import os
import sys

import healpy as hp
from mpi4py import MPI
import numpy as np
import matplotlib.pyplot as plt


comm = MPI.COMM_WORLD
ntask = comm.size
rank = comm.rank

if rank == 0:
    print(f"Running with {ntask} MPI tasks")
prefix = f"{rank:04} : "

alt_bands = {
    "f020" : "f020",
    "f026" : "f030",
    "f039" : "f040",
    "f085" : "f085",
    "f092" : "f090",
    "f095" : "f095",
    "f145" : "f145",
    "f155" : "f155",
    "f148" : "f150",
    "f227" : "f220",
    "f286" : "f280",
}

indir = "/global/cfs/cdirs/cmbs4/chile_optimization/simulations/phase1/noise_10_years"

nmc = 4
# fsky = 0.03
fsky = 0.01

nyear_cov = 10  # ab initio survey length
nyear_sim = 10  # performance-based simulation survey length
nyear_target = 10  # scale both to this survey length

ijob = -1
cldir = "noise_cl_phase1"
os.makedirs(cldir, exist_ok=True)
for config in "SAT",:
    for first_file in sorted(glob.glob(f"{indir}/phase1_noise_*_{config}_mc_0000.fits")):
        ijob += 1
        if ijob % ntask != rank:
            continue
        print(prefix + f"Processing {first_file}")

        band = os.path.basename(first_file).split("_")[2]
        alt_band = alt_bands[band]
        fname_cov = glob.glob(f"/global/cfs/cdirs/cmbs4/chile_optimization/simulations/phase1/noise_depth/sat_{alt_band}_*years_cov.fits")[0]
        # print(f"Loading {fname_cov}")
        cov = hp.read_map(fname_cov, [0, 3, 5])
        cov *= nyear_cov / nyear_target
        nside_in = hp.get_nside(cov)
        nside_out = 512
        cov = hp.ud_grade(cov, nside_out) * (nside_out / nside_in)**2

        basename = os.path.basename(first_file)
        outfile = os.path.join(cldir, basename.replace("_0000.fits", "_cl.fits"))
        outfile_white = os.path.join(cldir, basename.replace("_0000.fits", "_cl_white.fits"))
        if not os.path.isfile(outfile):
            bad = cov[0] == 0
            noisevar = cov[1] + cov[2]
            noisevar[bad] = 1e30
            limit = np.sort(noisevar)[int(fsky * noisevar.size)]
            mask = 1 - hp.smoothing(noisevar < limit, fwhm=np.radians(3), lmax=256, iter=0)
            limit = np.sort(mask)[int(fsky * mask.size)]
            mask = mask < limit
            maskfile = os.path.join(cldir, basename.replace("_0000.fits", "_mask.fits"))
            if maskfile == first_file:
                raise RuntimeError("Bad output file")
            hp.write_map(maskfile, mask, dtype=bool, overwrite=True)
            print(prefix + f"Wrote {maskfile}")

            mask_cov = hp.ud_grade(mask * 1.0, hp.get_nside(cov)) > 0.5
            cov[:, np.logical_not(mask_cov)] = 0

            clmean = None
            clmean_white = None
            for mc in range(nmc):
                np.random.seed(int(band[1:]) + 1000 * mc + 1000000 * ijob)
                white_noise = np.random.randn(cov.size).reshape(cov.shape) * np.sqrt(cov)
                fname_in = first_file.replace("_0000.fits", f"_{mc:04}.fits")
                fname_out = fname_in.replace(indir, cldir)
                fname_white = fname_out.replace(".fits", "_white.fits")
                if os.path.isfile(fname_out):
                    print(prefix + f"Reading {fname_out}")
                    cl = hp.read_cl(fname_out)
                    cl_white = hp.read_cl(fname_white)
                else:
                    print(prefix + f"    Loading {fname_in}")
                    m = hp.read_map(fname_in, None, dtype=np.float64)
                    m[:, np.logical_not(mask)] = 0
                    m *= np.sqrt(nyear_sim / nyear_target)
                    print(prefix + f"    Computing {fname_out}")
                    cl = hp.anafast(m, lmax=3 * hp.get_nside(m), iter=0) / fsky
                    cl_white = hp.anafast(white_noise, lmax=3 * hp.get_nside(white_noise), iter=0) / fsky
                    print(prefix + f"Writing {fname_out}")
                    hp.write_cl(fname_out, cl)
                    hp.write_cl(fname_white, cl_white)
                if clmean is None:
                    clmean = cl
                    clmean_white = cl_white
                else:
                    clmean += cl
                    clmean_white += cl_white
            clmean /= nmc
            clmean_white /= nmc
            print(prefix + f"    Writing {outfile}")
            hp.write_cl(outfile, clmean)
            hp.write_cl(outfile_white, clmean_white)

        clmean = hp.read_cl(outfile)
        clmean_white = hp.read_cl(outfile_white)

        # Plot the BB spectrum
        def noise_depth(level):
            # Convert noise level in (K.rad)^2 to uK.arcmin
            return np.sqrt(level) * 1e6 * 180 / np.pi * 60
        clmean = clmean[2]
        clmean_white = clmean_white[2]
        level = np.mean(clmean[900:1000])
        depth = noise_depth(level)
        level_white = np.mean(clmean_white[500:1000])
        depth_white = noise_depth(level_white)
        nrow, ncol = 1, 1
        fig = plt.figure(figsize=[6 * ncol, 4 * nrow])
        ax = fig.add_subplot(nrow, ncol, 1)
        scale = 1e12
        ax.plot(clmean * scale, "-", color="tab:blue", label="Sim")
        ax.plot(clmean_white * scale, "-", color="tab:orange", label="White noise")
        ax.axhline(level * scale, ls="--", color="tab:blue", label=f"{depth:.3} uK.arcmin")
        ax.axhline(level_white * scale, ls="--", color="tab:orange", label=f"{depth_white:.3} uK.arcmin")
        ax.legend(loc="best")
        ax.set_xlabel(r"Multipole, $\ell$")
        ax.set_ylabel(r"C$_\ell^\mathrm{BB}$ [$\mu$K$^2$]")
        ratio = depth / depth_white
        ax.set_title(f"{config} {band} BB, sim/depth map = {ratio:.3f} = SQRT(1 / {ratio**(-2):.3f}), fsky={fsky}")
        fname_plot = f"{cldir}/cl_comparison_{config}_{band}.png"
        fig.savefig(fname_plot)
        print(f"Wrote {fname_plot}")

comm.Barrier()
if rank == 0:
    print("All done!")
comm.Barrier()
