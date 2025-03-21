# 2025-03-20: Original version by Shamik Ghosh

import glob
import os
import sys

import astropy.units as u
import healpy as hp
from mpi4py import MPI
import numpy as np
import os

comm = MPI.COMM_WORLD
ntask = comm.size
rank = comm.rank

if rank == 0:
    print(f"Running with {ntask} MPI tasks")
prefix = f"{rank:04} : "

npipe_noise_dir  = '/global/cfs/cdirs/cmb/data/planck2020/npipe/npipe6v20_sim/residual/'
beam_dir         = '/global/cfs/cdirs/cmb/data/planck2020/npipe/npipe6v20/quickpol/'
fgdir            = "/global/cfs/cdirs/cmb/gsharing/panexp_v1_planck/galactic_foregrounds_mediumcomplexity/"
cmbdir           = "/global/cfs/cdirs/cmb/data/generic/cmb/ffp10/mc/scalar/"
tensordir        = "/global/cfs/cdirs/cmb/data/generic/cmb/ffp10/mc/tensor/"
outputdir        = "/pscratch/sd/s/shamikg/chile_optimization_sims/planck/"

os.makedirs(outputdir, exist_ok=True)
os.makedirs(f"{outputdir}cmb/", exist_ok=True)
os.makedirs(f"{outputdir}noise/", exist_ok=True)
os.makedirs(f"{outputdir}total/", exist_ok=True)

npipe_noise_ini = 200

freqs = [
            "030",
            "044",
            "070",
            "100",
            "143",
            "217",
            "353",
            "545",
            "857",
]

rot = hp.Rotator(coord=('G', 'C'))

ijob = -1
for band in freqs:
    
    print(prefix + f"band = {band} GHz")

    # Get the foregrounds

    nside = 2048
    lmin  = 0
    lmax  = 2 * nside
    
    fname_Bl = f"{beam_dir}/Bl_TEB_npipe6v20_{band}GHzx{band}GHz.fits"
    if int(band) < 100:
        fname_fg = f"{fgdir}" \
            f"sobs_mbs-s0017-20250208_LFI_mission_{band}_galactic_foregrounds_mediumcomplexity_healpix.fits"
    else:
        fname_fg = f"{fgdir}" \
            f"sobs_mbs-s0017-20250208_HFI_mission_{band}_galactic_foregrounds_mediumcomplexity_healpix.fits"

    fg = None
    Bl = None
    for mc in range(100):

        ijob += 1
        if ijob % ntask != rank:
            continue

        # Check if the last maps for this MC were already written
        if int(band) < 100:
            fname_total = f"{outputdir}total/planck_total_{band}_mc_{mc:04}.fits"
            if os.path.isfile(fname_total):
                print(prefix + f"    {fname_total} exists, skipping mc = {mc}")
                continue
            else:
                print(prefix + f"    {fname_total} does not exist.")
        elif int(band) >= 100:
            fname_total5 = f"{outputdir}total/planck_total_{band}_mc_{mc:04}.fits"
            if os.path.isfile(fname_total5):
                print(prefix + f"    {fname_total5} exists, skipping mc = {mc}")
                continue

        # Proceed with the simulation

        if fg is None:
            # First time foregrounds are needed
            print(prefix + f"    Reading {fname_fg}")
            fg = hp.read_map(fname_fg, field=None) * 1e-6  # to K_CMB
            print(prefix + f"    Expanding foreground in alm")
            fg_alms = hp.map2alm(fg, lmax=lmax, use_pixel_weights=True)
            
            
        # Get the beam
        if Bl is None:
            print(prefix + f"    Reading {fname_Bl}")
            Bl = hp.read_cl(fname_Bl)
            

        # Process all survey lengths for this band and this MC

        # Load the CMB realization

        print(prefix + f"    MC = {mc}")
        fname_cmb = f"{cmbdir}/ffp10noaber_lensed_scl_cmb_000_alm_mc_{mc:04}.fits"
        print(prefix + f"        Reading {fname_cmb}")
        alms = hp.read_alm(fname_cmb, hdu=(1, 2, 3))
        lmax_in = hp.Alm.getlmax(alms[0].size)
        if lmax != lmax_in:
            alms = hp.resize_alm(alms, lmax_in, lmax_in, lmax, lmax)

        # Add primordial B-modes to even-numbered MC

        # if mc % 2 == 0:
        #     fname_tensor = f"{tensordir}/ffp10_ten_cmb_000_alm_mc_{mc:04}.fits"
        #     print(prefix + f"        Reading {fname_tensor}")
        #     tensor_alms = hp.read_alm(fname_tensor, hdu=(1, 2, 3))
        #     lmax_in = hp.Alm.getlmax(tensor_alms[0].size)
        #     if lmax != lmax_in:
        #         tensor_alms = hp.resize_alm(tensor_alms, lmax_in, lmax_in, lmax, lmax)
        #     r_in = 0.01  # FFP10 simulations have r=0.01
        #     r_out = 0.003  # We want r = 0.003
        #     scale = np.sqrt(r_out / r_in)  # Power spectrum ratios in map domain
        #     for i in range(3):
        #         alms[i] += scale * tensor_alms[i]
        #     r = r_out
        # else:
        fname_tensor = None
        r = 0

        for i in range(3):
            hp.almxfl(alms[i], Bl[i], inplace=True)

        # Make a pure CMB map for reference

        print(prefix + "        Synthesizing CMB map")
        cmb = hp.alm2map(alms, nside)  # K_CMB
        fname_cmb = f"{outputdir}cmb/planck_cmb_{band}_mc_{mc:04}.fits"
        print(prefix + f"        Writing {fname_cmb}")
        os.makedirs(os.path.dirname(fname_cmb), exist_ok=True)
        hp.write_map(
            fname_cmb,
            cmb.astype(np.float32),
            dtype=np.float32,
            coord="C",
            column_units="K_CMB",
            extra_header=[
                ("BEAM", f"NPIPE TEB {band}GHzx{band}GHz"),
                ("CMB", fname_cmb),
                ("TENSOR", fname_tensor),
                ("R", r, "tensor-scalar ratio"),
            ],
            overwrite=True,
        )

        # Combine the foregrounds with the CMB and translate into a map
        # UPDATE 2024-11-20: Beam smoothing and pixel window smoothing are now not
        #                    performed on the total signal alms.
        for i in range(3):
            alms[i] += fg_alms[i]
        print(prefix + "        Synthesizing sky map")
        sky = hp.alm2map(alms, nside)  # K_CMB


        # Read and rotate noise maps
        noise_mc = npipe_noise_ini+mc
        fname_npipe = f"{npipe_noise_dir}{noise_mc:04}residual_npipe6v20_{band}_{noise_mc:04}.fits"
        print(prefix + f"        Reading {fname_npipe}")
        noise_npipe = hp.read_map(fname_npipe, field=None)
        noise_npipe[noise_npipe == hp.UNSEEN] = 0.
        noise_npipe = rot.rotate_map_pixel(noise_npipe)
        
        fname_noise = f"{outputdir}noise/planck_noise_{band}_mc_{mc:04}.fits"
        print(prefix + f"        Writing {fname_noise}")
        os.makedirs(os.path.dirname(fname_noise), exist_ok=True)

        args = {
                    "dtype" : np.float32,
                    "coord" : "C",
                    "column_units" : "K_CMB",
                    "extra_header" : [
                        ("NOISE", f"Planck NPIPE {band}"),
                        ("NPIPE MC", noise_mc),
                        ("NPIPE FILE", os.path.basename(fname_npipe)),
                    ],
                    "overwrite" : True,
                }
        hp.write_map(fname_noise, noise_npipe.astype(np.float32), **args)
        
        
        # coadd all maps and save to file
        total = sky + noise_npipe
        fname_total = f"{outputdir}total/planck_total_{band}_mc_{mc:04}.fits"
        os.makedirs(os.path.dirname(fname_total), exist_ok=True)
        args = {
            "dtype" : np.float32,
            "coord" : "C",
            "column_units" : "K_CMB",
            "extra_header" : [
                ("BEAM", f"NPIPE TEB {band}GHzx{band}GHz"),
                ("FOREGRND", os.path.basename(fname_fg)),
                ("CMB", os.path.basename(fname_cmb)),
                ("TENSOR", os.path.basename(fname_tensor)),
                ("R", r, "tensor-scalar ratio"),
                ("LMIN_CMB", lmin, "High-pass cut-off"),
                ("NOISE", os.path.basename(fname_noise)),
            ],
            "overwrite" : True,
        }
        print(prefix + f"        Writing {fname_total}")
        hp.write_map(fname_total, total.astype(np.float32), **args)

comm.Barrier()
if rank == 0:
    print("All done!")
comm.Barrier()



