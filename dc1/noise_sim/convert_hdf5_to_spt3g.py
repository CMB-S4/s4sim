#!/usr/bin/env python3

"""
This script loads HDF5 format data and exports to SPT3G format data.
"""

import os
import sys
import shutil
import re
import glob

import datetime

import argparse

import numpy as np

from astropy import units as u

import toast
import toast.ops

from toast.timing import gather_timers, dump

from toast import spt3g as t3g
from spt3g import core as c3g

from toast.observation import default_values as defaults


def parse_arguments():
    """
    Defines and parses the arguments for the script.
    """
    parser = argparse.ArgumentParser(
        description="Convert CMB-S4 simulated HDF5 data to SPT3G format"
    )

    parser.add_argument(
        "--obs",
        type=str,
        required=False,
        nargs="+",
        help="One or more observation directories.",
    )

    # The operators we want to configure from the command line or a parameter file.
    operators = list()

    # Parse all of the operator configuration
    config, args, jobargs = toast.parse_config(parser, operators=operators)

    return config, args, jobargs


def main():
    env = toast.utils.Environment.get()
    log = toast.utils.Logger.get()
    env.enable_function_timers()
    global_timer = toast.timing.GlobalTimers.get()
    global_timer.start("convert_hdf5_to_spt3g (total)")

    config, args, jobargs = parse_arguments()

    # Use a group size of one
    comm = toast.Comm(groupsize=jobargs.group_size)

    # Create the (initially empty) data
    data = toast.Data(comm=comm)

    # Set up temporary symlinks and output
    nowstr = datetime.datetime.now(datetime.timezone.utc).isoformat()
    input = f"hdf5_{nowstr}"
    output = f"spt3g_{nowstr}"

    if comm.world_rank == 0:
        os.makedirs(input)
        os.makedirs(output)

    obsloc = dict()

    for odir in args.obs:
        obsdir = os.path.basename(os.path.normpath(odir))
        # Look for the data file
        obsdata = glob.glob(os.path.join(odir, "RISING*.h5"))
        if len(obsdata) == 0:
            obsdata = glob.glob(os.path.join(odir, "SETTING*.h5"))
        if len(obsdata) == 0:
            obsdata = glob.glob(os.path.join(odir, "POLE*.h5"))
        if len(obsdata) != 1:
            raise RuntimeError(f"Did not find exactly one data file for {odir}")
        obsdata = os.path.abspath(obsdata[0])
        obsdatafile = os.path.basename(obsdata)
        obsname = re.sub(r"_\d+\.h5", "", obsdatafile)
        obsloc[obsname] = obsdir
        log.info_rank(
            f"Found observation data {obsdatafile}",
            comm=comm.comm_world
        )
        if comm.world_rank == 0:
            os.symlink(obsdata, os.path.join(input, obsdatafile))

    # Load the data.
    loader = toast.ops.LoadHDF5(
        volume=input,
        process_rows=1,
    )
    log.info_rank(
        f"Loading HDF5 data from {loader.volume}", comm=comm.comm_world
    )
    loader.apply(data)

    # Build up the lists of objects to export from the first observation

    noise_models = list()
    meta_arrays = list()
    shared = list()
    detdata = list()
    intervals = list()

    msg = "Exporting observation fields:"

    ob = data.obs[0]
    for k, v in ob.shared.items():
        g3name = f"shared_{k}"
        if re.match(r".*boresight.*", k) is not None:
            # These are quaternions
            msg += f"\n  (shared):    {k} (quaternions)"
            shared.append((k, g3name, c3g.G3VectorQuat))
        elif k == defaults.times:
            # Timestamps are handled separately
            continue
        else:
            msg += f"\n  (shared):    {k}"
            shared.append((k, g3name, None))
    for k, v in ob.detdata.items():
        msg += f"\n  (detdata):   {k}"
        if k == "signal":
            # We are going to convert this below.  Append
            # the new name to the list.
            detdata.append(("signal_f32", k, None))
        else:
            detdata.append((k, k, None))
    for k, v in ob.intervals.items():
        msg += f"\n  (intervals): {k}"
        intervals.append((k, k))
    for k, v in ob.items():
        if isinstance(v, toast.noise.Noise):
            msg += f"\n  (noise):     {k}"
            noise_models.append((k, k))
        elif isinstance(v, np.ndarray) and len(v.shape) > 0:
            if isinstance(v, u.Quantity):
                raise NotImplementedError("Writing array quantities not yet supported")
            msg += f"\n  (meta arr):  {k}"
            meta_arrays.append((k, k))
        else:
            msg += f"\n  (meta):      {k}"

    log.info_rank(msg, comm=comm.comm_world)

    # Convert the 8-byte detector data to 4-byte.

    log.info_rank("Converting signal to 4-byte floats", comm=comm.comm_world)
    for ob in data.obs:
        ob.detdata.create("signal_f32", dtype=np.float32)
        for d in ob.local_detectors:
            ob.detdata["signal_f32"][d, :] = ob.detdata["signal"][d, :]
        del ob.detdata["signal"]

    # Export the data

    meta_exporter = t3g.export_obs_meta(
        noise_models=noise_models,
        meta_arrays=meta_arrays,
    )
    data_exporter = t3g.export_obs_data(
        shared_names=shared,
        det_names=detdata,
        interval_names=intervals,
        compress=True,
    )
    exporter = t3g.export_obs(
        meta_export=meta_exporter,
        data_export=data_exporter,
        export_rank=0,
    )

    save3g = toast.ops.SaveSpt3g(
        directory=output,
        framefile_mb=1000,
        obs_export=exporter,
        purge=True,
    )

    log.info_rank(
        f"Exporting SPT3G data to {save3g.directory}", comm=comm.comm_world
    )
    save3g.apply(data)

    if comm.world_rank == 0:
        # Move spt3g outputs to the observation directories
        for out3g in obsloc.keys():
            os.rename(os.path.join(output, out3g), os.path.join(obsloc[out3g], "spt3g"))
        # Remove the input symlinks
        shutil.rmtree(input)

    # Dump all the timing information to the output dir

    global_timer.stop("convert_hdf5_to_spt3g (total)")
    alltimers = gather_timers(comm=comm.comm_world)
    if comm.world_rank == 0:
        out = os.path.join(output, "timing")
        dump(alltimers, out)


if __name__ == "__main__":
    world, procs, rank = toast.mpi.get_world()
    with toast.mpi.exception_guard(comm=world):
        main()
