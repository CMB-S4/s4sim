#!/usr/bin/env python3

"""
This script loads v0 uncompressed TOAST observations and writes v1 compressed files.
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

from toast.timing import gather_timers, dump, Timer

from toast.observation import default_values as defaults


def parse_arguments():
    """
    Defines and parses the arguments for the script.
    """
    parser = argparse.ArgumentParser(
        description="Compress CMB-S4 simulation data"
    )

    parser.add_argument(
        "--obs",
        type=str,
        required=False,
        nargs="+",
        help="One or more observation files",
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
    global_timer.start("compress HDF5 (total)")

    config, args, jobargs = parse_arguments()

    # Default group size
    comm = toast.Comm()

    # Process each observation
    for obs_path in args.obs:
        obs_dir = os.path.dirname(obs_path)
        file_root = os.path.splitext(obs_path)[0]
        if comm.world_rank == 0:
            print(f"Working on {obs_path}:")
        backup = f"{file_root}_uncompressed.h5"
        timer = Timer()
        timer.start()
        obs = toast.io.load_hdf5(
            obs_path,
            comm,
            process_rows=comm.group_size,
            meta=None,
            detdata=None,
            shared=None,
            intervals=None,
            detectors=None,
            force_serial=False,
        )

        if comm.comm_world is not None:
            comm.comm_world.barrier()
        timer.stop()
        if comm.world_rank == 0:
            print(f"  Load {obs_path} in {timer.seconds()} s", flush=True)

        if comm.world_rank == 0:
            os.rename(obs_path, backup)
        
        if comm.comm_world is not None:
            comm.comm_world.barrier()
        
        timer.start()
        obf = toast.io.save_hdf5(
            obs,
            obs_dir,
            meta=None,
            detdata=[
                (defaults.det_data, {"type": "flac"}),
                (defaults.det_flags, {"type": "gzip"}),
            ],
            shared=None,
            intervals=None,
            config=None,
            times=defaults.times,
            force_serial=False,
        )
        if comm.comm_world is not None:
            comm.comm_world.barrier()
        timer.stop()
        if comm.world_rank == 0:
            print(f"  Save {obs_path} in {timer.seconds()} s", flush=True)

        if obf != obs_path:
            msg = f"Generated HDF5 ({obf}) does not match original "
            msg += f"file name ({obs_path})"
            raise RuntimeError(msg)

        timer.start()
        compare = toast.io.load_hdf5(
            obs_path,
            comm,
            process_rows=comm.group_size,
        )
        if comm.comm_world is not None:
            comm.comm_world.barrier()
        timer.stop()
        if comm.world_rank == 0:
            print(
                f"  Re-load {obs_path} for verification in {timer.seconds()} s",
                flush=True
            )

        if compare != obs:
            msg = f"Observation HDF5 verify failed:\n"
            msg += f"Input = {obs}\n"
            msg += f"Loaded = {compare}"
            log.error(msg)
            raise RuntimeError(msg)
        elif comm.world_rank == 0:
            print(f"  Verification PASS", flush=True)

    # Dump all the timing information to the output dir

    global_timer.stop("compress HDF5 (total)")
    alltimers = gather_timers(comm=comm.comm_world)
    if comm.world_rank == 0:
        out = os.path.join(".", "timing")
        dump(alltimers, out)


if __name__ == "__main__":
    world, procs, rank = toast.mpi.get_world()
    with toast.mpi.exception_guard(comm=world):
        main()
