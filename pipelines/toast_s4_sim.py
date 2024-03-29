#!/usr/bin/env python

# Copyright (c) 2020-2020 CMB-S4 Collaboration..
# Full license can be found in the top level "LICENSE" file.


import os

if "TOAST_STARTUP_DELAY" in os.environ:
    import numpy as np
    import time

    delay = np.float(os.environ["TOAST_STARTUP_DELAY"])
    wait = np.random.rand() * delay
    # print('Sleeping for {} seconds before importing TOAST'.format(wait),
    #      flush=True)
    time.sleep(wait)


# TOAST must be imported before numpy to ensure the right MKL is used
import toast

# import so3g

import copy
from datetime import datetime
import gc
import os
import pickle
import re
import sys
import traceback
from time import time

import argparse
import dateutil.parser

from toast.mpi import get_world, Comm
from toast.utils import Logger, Environment, memreport

from toast.timing import function_timer, GlobalTimers, Timer, gather_timers
from toast.timing import dump as dump_timing

import toast.pipeline_tools as toast_tools

import s4sim.pipeline_tools as s4_tools

import numpy as np

import s4sim.hardware

import warnings

warnings.filterwarnings("ignore")
# warnings.filterwarnings('error')
# warnings.simplefilter('ignore', ImportWarning)
# warnings.simplefilter('ignore', ResourceWarning)
# warnings.simplefilter('ignore', DeprecationWarning)
# warnings.filterwarnings("ignore", message="numpy.dtype size changed")
# warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

XAXIS, YAXIS, ZAXIS = np.eye(3)


def parse_arguments(comm):
    timer = Timer()
    timer.start()
    log = Logger.get()

    parser = argparse.ArgumentParser(
        description="Simulate ground-based boresight pointing.  Simulate "
        "atmosphere and make maps for some number of noise Monte Carlos.",
        fromfile_prefix_chars="@",
    )

    toast_tools.add_dist_args(parser)
    toast_tools.add_todground_args(parser)
    toast_tools.add_pointing_args(parser)
    toast_tools.add_polyfilter_args(parser)
    toast_tools.add_polyfilter2D_args(parser)
    toast_tools.add_common_mode_filter_args(parser)
    toast_tools.add_groundfilter_args(parser)
    toast_tools.add_atmosphere_args(parser)
    toast_tools.add_noise_args(parser)
    toast_tools.add_gainscrambler_args(parser)
    toast_tools.add_madam_args(parser)
    toast_tools.add_filterbin_args(parser)
    toast_tools.add_sky_map_args(parser)
    toast_tools.add_sss_args(parser)
    toast_tools.add_tidas_args(parser)
    toast_tools.add_mc_args(parser)
    s4_tools.add_hw_args(parser)
    s4_tools.add_s4_noise_args(parser)
    s4_tools.add_pysm_args(parser)
    toast_tools.add_debug_args(parser)

    parser.add_argument(
        "--no-maps",
        required=False,
        default=False,
        action="store_true",
        help="Disable all mapmaking.",
    )

    parser.add_argument(
        "--skip-madam",
        required=False,
        default=False,
        action="store_true",
        help="Skip the first Madam call.",
    )

    parser.add_argument(
        "--pairdiff",
        required=False,
        default=False,
        action="store_true",
        help="Pair-difference TOD and pointing.",
    )

    parser.add_argument(
        "--outdir", required=False, default="out", help="Output directory"
    )

    try:
        args = parser.parse_args()
    except SystemExit as e:
        sys.exit()

    if len(args.bands.split(",")) != 1:
        # Multi frequency run.  We don't support multiple copies of
        # scanned signal.
        if args.input_map:
            raise RuntimeError(
                "Multiple frequencies are not supported when scanning from a map"
            )

    if args.simulate_atmosphere and args.weather is None:
        raise RuntimeError("Cannot simulate atmosphere without a TOAST weather file")

    if comm.world_rank == 0:
        log.info("\n")
        log.info("All parameters:")
        for ag in vars(args):
            log.info("{} = {}".format(ag, getattr(args, ag)))
        log.info("\n")

    if args.group_size:
        comm = Comm(groupsize=args.group_size)

    if comm.world_rank == 0:
        os.makedirs(args.outdir, exist_ok=True)
        timer.report_clear("Parse arguments")

    return args, comm


def setup_output(args, comm, mc):
    outpath = "{}/{:08}".format(args.outdir, mc)
    if comm.world_rank == 0:
        os.makedirs(outpath, exist_ok=True)
    return outpath


def outputs_exist(args, comm, outpath):
    """ Returns True if all of the requested outputs already exist
    """
    there = True
    if comm.world_rank == 0:
        if not args.skip_madam:
            if there and args.write_binmap:
                fname = os.path.join(
                    outpath, args.mapmaker_prefix + "_telescope_all_time_all_bmap.fits"
                )
                there = os.path.isfile(fname)
                if there:
                    print(f"{fname} exists", flush=True)
                else:
                    print(f"{fname} does not exist", flush=True)
            if there and args.destripe:
                fname = os.path.join(
                    outpath, args.mapmaker_prefix + "_telescope_all_time_all_map.fits"
                )
                there = os.path.isfile(fname)
                if there:
                    print(f"{fname} exists", flush=True)
                else:
                    print(f"{fname} does not exist", flush=True)
        if there and (
                args.apply_polyfilter
                or args.apply_polyfilter2D
                or args.apply_common_mode_filter
                or args.apply_groundfilter
        ):
            fname = os.path.join(
                outpath,
                args.mapmaker_prefix + "_filtered" + "_telescope_all_time_all_bmap.fits",
            )
            there = os.path.isfile(fname)
            if there:
                print(f"{fname} exists", flush=True)
            else:
                print(f"{fname} does not exist", flush=True)
        if there and (args.filterbin_ground_order or args.filterbin_poly_order):
            fname = os.path.join(
                outpath,
                args.filterbin_prefix + "_telescope_all_time_all_filtered.fits",
            )
            there = os.path.isfile(fname) or os.path.isfile(fname + ".gz")
            if there:
                print(f"{fname} exists", flush=True)
            else:
                print(f"{fname} does not exist", flush=True)
    there = comm.comm_world.bcast(there)
    return there


def pairdiff(data, args, comm, name, do_pointing):
    if not args.pairdiff:
        return
    t1 = time()
    if comm.comm_world.rank == 0:
        print("Pair differencing data", flush=True)

    for obs in data.obs:
        tod = obs["tod"]
        for det in tod.local_dets:
            signal = tod.local_signal(det, name)
            if det.endswith("A"):
                pairdet = det[:-1] + "B"
                if pairdet not in tod.local_dets:
                    raise RuntimeError(
                        "Detector pair not available ({}, {})".format(det, pairdet)
                    )
            else:
                continue
            # signal
            pairsignal = tod.local_signal(pairdet, name)
            signal[:], pairsignal[:] = [
                0.5 * (signal + pairsignal), 0.5 * (signal - pairsignal)
            ]
            if do_pointing:
                # flags
                flags = tod.local_flags(det)
                pairflags = tod.local_flags(pairdet)
                flags |= pairflags
                pairflags[:] = flags
                # pointing weights
                weights = tod.cache.reference("weights_" + det)
                pairweights = tod.cache.reference("weights_" + pairdet)
                weights[:], pairweights[:] = [
                    0.5 * (weights + pairweights), 0.5 * (weights - pairweights)
                ]

    if comm.comm_world.rank == 0:
        print("Pair differenced in {:.1f} s".format(time() - t1), flush=True)
    return

def main():
    log = Logger.get()
    gt = GlobalTimers.get()
    gt.start("toast_s4_sim (total)")
    timer0 = Timer()
    timer0.start()

    mpiworld, procs, rank, comm = toast_tools.get_comm()

    memreport("at the beginning of the pipeline", comm.comm_world)

    args, comm = parse_arguments(comm)

    # Initialize madam parameters

    madampars = toast_tools.setup_madam(args)

    # Load and broadcast the schedule file

    schedules = toast_tools.load_schedule(args, comm)

    # Load the weather and append to schedules

    toast_tools.load_weather(args, comm, schedules)

    # load or simulate the focalplane

    detweights = s4_tools.load_focalplanes(args, comm, schedules)

    # Create the TOAST data object to match the schedule.  This will
    # include simulating the boresight pointing.

    data, telescope_data = s4_tools.create_observations(args, comm, schedules)

    memreport("after creating observations", comm.comm_world)

    # Optionally rewrite the noise PSD:s in each observation to include
    # elevation-dependence
    s4_tools.get_elevation_noise(args, comm, data)

    totalname = "total"

    # Split the communicator for day and season mapmaking

    time_comms = toast_tools.get_time_communicators(args, comm, data)

    # Expand boresight quaternions into detector pointing weights and
    # pixel numbers

    toast_tools.expand_pointing(args, comm, data)

    # Only purge the pointing if we are NOT going to export the
    # data to a TIDAS volume
    if args.tidas is None:
        for ob in data.obs:
            tod = ob["tod"]
            try:
                tod.free_radec_quats()
            except AttributeError:
                # These TOD objects do not have RA/Dec quaternions
                pass

    memreport("after pointing", comm.comm_world)

    # Prepare auxiliary information for distributed map objects

    memreport("after submaps", comm.comm_world)

    # Set up objects to take copies of the TOD at appropriate times

    if args.pysm_model:
        if schedules is not None:
            focalplanes = [s.telescope.focalplane.detector_data for s in schedules]
        else:
            focalplanes = [telescope.focalplane.detector_data]
        signalname = s4_tools.simulate_sky_signal(args, comm, data, focalplanes)
    else:
        signalname = toast_tools.scan_sky_signal(args, comm, data)

    memreport("after PySM", comm.comm_world)

    # Loop over Monte Carlos

    firstmc = int(args.MC_start)
    nmc = int(args.MC_count)

    for mc in range(firstmc, firstmc + nmc):

        if comm.world_rank == 0:
            log.info("Processing MC = {}".format(mc))

        # Uncomment to run with new TOAST
        #toast_tools.draw_weather(args, comm, data, mc)

        outpath = setup_output(args, comm, mc)

        if outputs_exist(args, comm, outpath):
            if comm.world_rank == 0:
                log.info("Outputs already exist, skipping.")
            continue

        toast.tod.OpCacheClear(totalname).exec(data)

        toast_tools.simulate_atmosphere(args, comm, data, mc, totalname)

        s4_tools.scale_atmosphere_by_bandpass(args, comm, data, totalname, mc)

        memreport("after atmosphere", comm.comm_world)

        # update_atmospheric_noise_weights(args, comm, data, freq, mc)

        toast_tools.add_signal(
            args, comm, data, totalname, signalname, purge=(mc == firstmc + nmc - 1)
        )

        memreport("after adding sky", comm.comm_world)

        toast_tools.simulate_noise(args, comm, data, mc, totalname)

        memreport("after simulating noise", comm.comm_world)

        toast_tools.simulate_sss(args, comm, data, mc, totalname)

        memreport("after simulating SSS", comm.comm_world)

        toast_tools.scramble_gains(args, comm, data, mc, totalname)

        if mc == firstmc:
            # For the first realization and frequency, optionally
            # export the timestream data.
            toast_tools.output_tidas(args, comm, data, totalname)

            memreport("after export", comm.comm_world)

        if args.no_maps:
            continue

        # Bin and destripe maps

        pairdiff(data, args, comm, totalname, mc == firstmc)

        if not args.skip_madam:
            toast_tools.apply_madam(
                args,
                comm,
                data,
                madampars,
                outpath,
                detweights,
                totalname,
                time_comms=time_comms,
                telescope_data=telescope_data,
                first_call=(mc == firstmc),
            )
            memreport("after madam", comm.comm_world)

        if (
                args.filterbin_ground_order is not None
                or args.filterbin_poly_order is not None
        ):
            toast_tools.apply_filterbin(
                args,
                comm,
                data,
                outpath,
                totalname,
                time_comms=time_comms,
                telescope_data=telescope_data,
                first_call=(mc == firstmc),
            )

        if (
                args.apply_polyfilter
                or args.apply_polyfilter2D
                or args.apply_common_mode_filter
                or args.apply_groundfilter
        ):

            # Filter signal

            toast_tools.apply_common_mode_filter(args, comm, data, totalname)

            toast_tools.apply_polyfilter2D(args, comm, data, totalname)

            toast_tools.apply_polyfilter(args, comm, data, totalname)

            toast_tools.apply_groundfilter(args, comm, data, totalname)

            memreport("after filter", comm.comm_world)

            # Bin maps

            toast_tools.apply_madam(
                args,
                comm,
                data,
                madampars,
                outpath,
                detweights,
                totalname,
                time_comms=time_comms,
                telescope_data=telescope_data,
                first_call=(args.skip_madam and mc == firstmc),
                extra_prefix="filtered",
                bin_only=(not args.skip_madam),
            )

            memreport("after filter & bin", comm.comm_world)

    if comm.comm_world is not None:
        comm.comm_world.barrier()

    memreport("at the end of the pipeline", comm.comm_world)

    gt.stop_all()
    if mpiworld is not None:
        mpiworld.barrier()
    timer = Timer()
    timer.start()
    alltimers = gather_timers(comm=mpiworld)
    if rank == 0:
        out = os.path.join(args.outdir, "timing")
        dump_timing(alltimers, out)
        timer.stop()
        timer.report("Gather and dump timing info")
    timer0.stop()
    if comm.world_rank == 0:
        timer0.report("toast_s4_sim.py pipeline")
    return


if __name__ == "__main__":

    try:
        main()
    except Exception as e:
        # We have an unhandled exception on at least one process.  Print a stack
        # trace for this process and then abort so that all processes terminate.
        mpiworld, procs, rank = get_world()
        exc_type, exc_value, exc_traceback = sys.exc_info()
        lines = traceback.format_exception(exc_type, exc_value, exc_traceback)
        lines = ["Proc {}: {}".format(rank, x) for x in lines]
        print("".join(lines), flush=True)
        if mpiworld is not None and procs > 1:
            mpiworld.Abort(6)
        else:
            raise
