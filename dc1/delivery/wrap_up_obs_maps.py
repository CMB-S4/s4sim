import h5py
import os
import sys

import healpy as hp
import numpy as np


# Filenames are based on
# https://docs.google.com/document/d/1VIQYsGMza9rOn3E0GY1hDpzLN4hcR7Fwm7tpXGUxSlY/edit?usp=sharing

telescopes_to_bands = {
    "chlat" : ["025", "040", "090", "150", "230", "280"],
}

# Simulation scripts used different names than the delivery

alternate_names = {
    "chlat" : "LAT0_CHLAT",
    "025" : "f030",
    "040" : "f040",
    "090" : "f090",
    "150" : "f150",
    "230" : "f220",
    "280" : "f280",
    "cmb" : "primary CMB",
    "cmb_lensing" : "lensing perturbation",
    "foreground" : "extragalactic + galactic foregrounds",
    "noise" : "atmosphere + noise",
}

fname_schedule = "../scan_strategy/chile_lat/schedules/chile_schedule_lat.pruned.txt"
schedule = []
with open(fname_schedule, "r") as file_in:
    for iline, line in enumerate(file_in):
        if iline < 3:
            continue
        parts = line.split()
        start_date = parts[0]
        start_time = parts[1]
        stop_date = parts[2]
        stop_time = parts[3]
        start_mjd = parts[4]
        stop_mjd = parts[5]
        start = f"{start_date} {start_time}"
        stop = f"{stop_date} {stop_time}"
        target = parts[7]
        scan = parts[21]
        subscan = parts[22]
        observation = f"{target}-{scan}-{subscan}"
        schedule.append((observation, start, stop, start_mjd, stop_mjd))

# Subsets

n_wafer_split = 1
i_wafer_split = 1
wafer_split_name = f"w{n_wafer_split:02}.{i_wafer_split:02}"
realization = 0
realization_name = f"r{realization:03}"
components = ["noise", "foreground", "cmb_lensing", "cmb"]
# four-digit mask identifying signal content.  The position of each
# digit matches an entry in `components`
#complements = ["0001", "0010", "0100", "1000", "1111"]
complements = ["0001", "0010", "0100", "1000"]

# Which data products to build

signal_products = [
    ("map01", "noise-weighted filter+bin IQU map"),
]

supporting_products = [
    ("mat01", "inverse white noise covariance matrix"),
]

rootdir = "/global/cfs/cdirs/cmbs4/dc/dc1/staging"
outdir = "/global/cfs/cdirs/cmbs4/dc/dc1/delivery"

def add_chunked_dataset(dset1, dset2):
    """ Add chunks from dset2 into dset1 """
    if dset1.shape != dset2.shape:
        raise RuntimeError("Mismatch in dataset dimensions")
    chunksize = dset2.chunks
    nchunk = dset2.id.get_num_chunks()
    for ichunk in range(nchunk):
        info = dset2.id.get_chunk_info(ichunk)
        offset = info.chunk_offset
        ind = tuple([slice(off, off + sz) for (off, sz) in zip(offset, chunksize)])
        dset1[ind] += dset2[ind]
    return

def subtract_chunked_dataset(dset1, dset2):
    """ Subtract chunks from dset2 from dset1 """
    if dset1.shape != dset2.shape:
        raise RuntimeError("Mismatch in dataset dimensions")
    chunksize = dset2.chunks
    nchunk = dset2.id.get_num_chunks()
    for ichunk in range(nchunk):
        info = dset2.id.get_chunk_info(ichunk)
        offset = info.chunk_offset
        ind = tuple([slice(off, off + sz) for (off, sz) in zip(offset, chunksize)])
        dset1[ind] -= dset2[ind]
    return

for telescope, bands in telescopes_to_bands.items():
    alt_telescope = alternate_names[telescope]
    for band in bands:
        if len(sys.argv) > 1 and band not in sys.argv[1:]:
            continue
        alt_band = alternate_names[band]
        dir_out = f"{outdir}/dc0/observation/{telescope}/{band}"
        os.makedirs(dir_out, exist_ok=True)
        for iobs, obs in enumerate(schedule):
            obs_id, start, stop, start_mjd, stop_mjd = obs
            obs_name = f"o{iobs:05}"
            
            # Supporting products do not reference signal component

            supporting_metadata = {
                "wafer_split" : n_wafer_split,
                "i_wafer_split" : i_wafer_split,
                "realization" : realization,
                "telescope" : telescope,
                "band" : band,
                "start" : start,
                "stop" : stop,
                "start_mjd" : start_mjd,
                "stop_mjd" : stop_mjd,
            }

            for product, description in supporting_products:
                fname_out = os.path.join(
                    dir_out,
                    f"dc0_{telescope}_{obs_name}_{band}_{product}.hdf5",
                )
                if os.path.isfile(fname_out):
                    print(f"{fname_out} exists, skipping ...", flush=True)
                    continue
                print(f"\nAssembling {fname_out} -- {description}", flush=True)
                if product == "mat01":
                    # Inverse white noise covariance in HDF5 format
                    fname_in = os.path.join(
                        rootdir,
                        f"noise_sim/outputs_rk/{alt_telescope}/{alt_band}",
                        f"{obs_id}",
                        f"mapmaker_{obs_id}_invcov.h5",
                    )
                    if not os.path.isfile(fname_in):
                        print(f"ERROR: missing input file: {fname_in}", flush=True)
                        continue
                    print(f"Loading {fname_in}", flush=True)
                    with h5py.File(fname_out, "w") as dest:
                        with h5py.File(fname_in, "r") as source:
                            source.copy(source["map"], dest)
                        for key, value in supporting_metadata.items():
                            dest.attrs[key] = value
                        dest.attrs["product"] = product
                        dest.attrs["units"] = "K_CMB^-2"
                    print(f"Done!", flush=True)
                else:
                    msg = f"Don't know how to assemble supporting product: {product}"
                    raise RuntimeError(msg)

            # Signal products come in many flavors

            for complement in complements:
                complement_name = f"c{complement}"
                for product, description in signal_products:
                    fname_out = os.path.join(
                        dir_out,
                        f"dc0_{telescope}_{obs_name}_{band}_{product}_{complement_name}.hdf5",
                    )
                    if os.path.isfile(fname_out):
                        print(f"{fname_out} exists, skipping ...", flush=True)
                        continue
                    print(f"\nAssembling {fname_out} -- {description}", flush=True)
                    try:
                        with h5py.File(fname_out, "w") as dest:
                            component_names = None
                            for i, component in enumerate(components):
                                if complement[i] == "0":
                                    continue
                                alt_component = alternate_names[component]
                                if product == "map01":
                                    # Noise-weighted filter-and-bin map
                                    fname_in = os.path.join(
                                        rootdir,
                                        f"{component}_sim/outputs_rk/{alt_telescope}/{alt_band}",
                                        f"{obs_id}",
                                        f"mapmaker_{obs_id}_noiseweighted_map.h5",
                                    )
                                else:
                                    msg = f"Don't know how to assemble signal product: {product}"
                                    raise RuntimeError(msg)
                                if not os.path.isfile(fname_in):
                                    msg = f"missing input file: {fname_in}"
                                    raise FileNotFoundError(msg)
                                print(f"Loading {fname_in}", flush=True)
                                with h5py.File(fname_in, "r") as source:
                                    if component_names is None:
                                        source.copy(source["map"], dest)
                                        component_names = f"{alt_component}"
                                    else:
                                        # Must loop over populated chunks,
                                        # otherwise the new file is dense
                                        add_chunked_dataset(dest["map"], source["map"])
                                        component_names += f", {alt_component}"
                                if component == "foreground":
                                    # Foreground maps had the CMB included
                                    # but we want it out
                                    fname_cmb = fname_in.replace("foreground", "cmb")
                                    if not os.path.isfile(fname_cmb):
                                        msg = f"missing input file: {fname_cmb}"
                                        raise FileNotFoundError(msg)
                                    print(f"Loading and subtracting {fname_cmb}", flush=True)
                                    with h5py.File(fname_cmb, "r") as source:
                                        subtract_chunked_dataset(dest["map"], source["map"])
                            for key, value in supporting_metadata.items():
                                dest.attrs[key] = value
                            dest.attrs["product"] = product
                            dest.attrs["units"] = "K_CMB^-1"
                            dest.attrs["complement"] = complement
                            dest.attrs["content"] = component_names
                        print(f"Done!", flush=True)
                    except Exception as e:
                        print(
                            f"ERROR: failed to assemble {fname_out}: '{e}'", flush=True
                        )
                        if os.path.isfile(fname_out):
                            print(f"Deleting failed output file: {fname_out}")
                            os.remove(fname_out)
