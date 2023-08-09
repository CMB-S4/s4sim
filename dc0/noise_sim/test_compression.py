import os
import sys

from toast.mpi import MPI, Comm
import toast
from toast.ops import LoadHDF5, SaveHDF5

"""

Execute the script with

OMP_NUM_THREADS=8 srun -N 1 -n 32 -c 8 --cpu-bind=cores python test_compression.py

"""

log = toast.utils.Logger.get()
timer = toast.timing.Timer()
timer.start()

wcomm, procs, rank = toast.get_world()
toast_comm = toast.Comm(world=wcomm)

indir = "outputs/LAT0_CHLAT/f090/RISING_SCAN_40-10-7"
data = toast.Data(comm=toast_comm)

nside = 4096
nest = False

loaded = False
# Precision=9 fails because 32bit FLAC cannot support it.
for precision in (2, 3, 4, 5, 6, 7, 8):
    outdir = f"testdata_precision{precision}"
    if os.path.isdir(outdir):
        log.info_rank(f"{outdir} exists, skipping...")
        continue
    if not loaded:
        log.info_rank(f"Loading data from {indir}", comm=wcomm)
        loader = LoadHDF5(
            volume=indir,
            detdata=["signal", "flags"],
            shared=["times", "flags", "boresight_radec"],
            pattern="^(RISING|SETTING).*h5",
        )
        loader.exec(data)
        log.info_rank(f"Loaded data from {indir} in", comm=wcomm, timer=timer)
        loaded = True
    log.info_rank(f"Writing data to {outdir}", comm=wcomm)
    writer = SaveHDF5(
        detdata_float32=True,
        volume=outdir,
        compress_detdata=True,
        compress_precision=precision,
    )
    writer.exec(data)
    log.info_rank(f"Wrote data to {outdir} in", comm=wcomm, timer=timer)

quat_pointing = toast.ops.PointingDetectorSimple()
pixel_pointing = toast.ops.PixelsHealpix(
    detector_pointing=quat_pointing, nside=nside, nest=nest
)
weights = toast.ops.StokesWeights(
    mode="IQU",
    detector_pointing=quat_pointing,
)
binner = toast.ops.BinMap(
    pixel_pointing=pixel_pointing,
    stokes_weights=weights,
)
mapmaker = toast.ops.MapMaker(
    binning=binner,
    write_hits=True,
    write_binmap=True,
    write_rcond=False,
    write_invcov=False,
    write_cov=False,
    write_map=False,
    write_hdf5=True,
    output_dir=indir,
)

for precision in (None, 2, 3, 4, 5, 6, 7, 8):
    if precision is not None:
        indir = f"testdata_precision{precision}"
    fname_binmap = f"{indir}/MapMaker_binmap.h5"
    if os.path.isfile(fname_binmap):
        log.info_rank(f"{fname_binmap} exists, skipping...", comm=wcomm)
        continue
    tempdata = toast.Data(comm=toast_comm)
    log.info_rank(f"Loading data from {indir}", comm=wcomm)
    loader = LoadHDF5(
        volume=indir,
        pattern="^(RISING|SETTING).*h5",
        force_serial=True,
    )
    loader.apply(tempdata)
    mapmaker.output_dir = indir
    mapmaker.apply(tempdata)
    del tempdata
