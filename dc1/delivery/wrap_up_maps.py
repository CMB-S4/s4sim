import os
import sys

import healpy as hp
import numpy as np


# Filenames are based on the https://docs.google.com/document/d/1VIQYsGMza9rOn3E0GY1hDpzLN4hcR7Fwm7tpXGUxSlY/edit?usp=sharing

telescopes = {
    "chlat" : ["027", "039", "093", "145", "225", "278"],
}

alternate_names = {
    "chlat" : "LAT0_CHLAT",
    "027" : "f030",
    "039" : "f040",
    "093" : "f090",
    "145" : "f150",
    "225" : "f220",
    "278" : "f280",
}

time_splits = [32, 16, 8, 4, 2, 1]
n_wafer_split = 1
i_wafer_split = 1
wafer_split_name = f"w{n_wafer_split:02}.{i_wafer_split:02}"
realization = 0
realization_name = f"r{realization:03}"

components = ["cmb", "noise", "foreground"]
complements = ["001", "010", "100", "111"]

signal_products = [
    # "map01",  # noise-weighted filter+bin iqu map
    "map02",  # filter+bin iqu map
]

supporting_products = [
    # "map00",  # hit map
    "map03",  # tp depth map
    # "mat01",  # inverse white noise covariance matrix
    # "mat02",  # white noise covariance matrix
]

rootdir = "/global/cfs/cdirs/cmbs4/dc/dc1/staging"
outdir = "/global/cfs/cdirs/cmbs4/dc/dc1/delivery"

for telescope, bands in telescopes.items():
    alt_telescope = alternate_names[telescope]
    for band in bands:
        alt_band = alternate_names[band]
        dir_out = f"{outdir}/dc0/mission/{telescope}/{band}"
        os.makedirs(dir_out, exist_ok=True)
        for n_time_split in time_splits:
            for i_time_split in range(1, n_time_split + 1):
                time_split_name = f"t{n_time_split:02}.{i_time_split:02}"

                # Supporting products do not reference signal component

                supporting_header = [
                    ("HIERARCH n_time_split", n_time_split, "Time-wise split level"),
                    ("HIERARCH i_time_split", i_time_split, "Time-wise split index"),
                    ("HIERARCH n_wafer_split", n_wafer_split, "Wafer split level"),
                    ("HIERARCH i_wafer_split", i_wafer_split, "Wafer split index"),
                    ("HIERARCH realization", realization, "Realization index"),
                    ("HIERARCH telescope", telescope, "Telescope"),
                    ("HIERARCH band", band, "Frequency band"),
                ]

                for product in supporting_products:
                    fname_out = os.path.join(
                        dir_out,
                        f"dc0_{telescope}_{band}_{product}_{wafer_split_name}_"
                        f"{time_split_name}_{realization_name}.fits",
                    )
                    if os.path.isfile(fname_out):
                        print(f"{fname_out} exists, skipping ...", flush=True)
                        continue
                    print(f"\nAssembling {fname_out}", flush=True)
                    if product == "map03":
                        # T+P depth map is derived from the 3x3 white noise covariance
                        fname_in = os.path.join(
                            rootdir,
                            f"noise_sim/outputs_rk/coadd/{alt_telescope}",
                            f"coadd_{alt_telescope}_{alt_band}_"
                            f"{i_time_split:03}of{n_time_split:03}_cov.fits",
                        )
                        print(f"Loading {fname_in}", flush=True)
                        ii, qq, uu = hp.read_map(fname_in, [0, 3, 5])
                        root_area = np.sqrt(hp.nside2pixarea(hp.get_nside(ii), degrees=True)) * 60
                        t_depth = np.sqrt(ii) * 1e6 * root_area
                        p_depth = np.sqrt(qq + uu) * 1e6 * root_area
                        total = np.vstack([t_depth, p_depth])
                        column_names = ["T_DEPTH", "P_DEPTH"]
                        column_units = "uK_arcmin"
                    elif product == "mat01":
                        # Inverse white noise covariance
                        fname_in = os.path.join(
                            rootdir,
                            f"noise_sim/outputs_rk/coadd/{alt_telescope}"
                            f"coadd_{alt_telescope}_{alt_band}_"
                            f"{i_time_split:03}of{n_time_split:03}_invcov.fits",
                        )
                        print(f"Loading {fname_in}", flush=True)
                        total = hp.read_map(fname_in, None)
                        column_names = ["II", "IQ", "IU", "QQ", "QU", "UU"]
                        column_units = "K^-2_CMB"
                    elif product == "mat02":
                        # White noise covariance
                        fname_in = os.path.join(
                            rootdir,
                            f"noise_sim/outputs_rk/coadd/{alt_telescope}",
                            f"coadd_{alt_telescope}_{alt_band}_"
                            f"{i_time_split:03}of{n_time_split:03}_cov.fits",
                        )
                        print(f"Loading {fname_in}", flush=True)
                        total = hp.read_map(fname_in, None)
                        column_names = ["II", "IQ", "IU", "QQ", "QU", "UU"]
                        column_units = "K^2_CMB"
                    else:
                        msg = f"Don't know how to assemble supporting product: {product}"
                        raise RuntimeError(msg)
                    print(f"Writing {fname_out}", flush=True)
                    hp.write_map(
                        fname_out,
                        total,
                        nest=False,
                        dtype=np.float32,
                        coord="C",
                        column_names=column_names,
                        column_units=column_units,
                        extra_header=supporting_header,
                    )

                # Signal products come in many flavors

                for complement in complements:
                    complement_name = f"c{complement}"
                    for product in signal_products:
                        fname_out = os.path.join(
                            dir_out,
                            f"dc0_{telescope}_{band}_{product}_{complement_name}_"
                            f"{wafer_split_name}_{time_split_name}_"
                            f"{realization_name}.fits",
                        )
                        if os.path.isfile(fname_out):
                            print(f"{fname_out} exists, skipping ...")
                            continue
                        print(f"\nAssembling {fname_out}", flush=True)
                        total = None
                        component_names = None
                        for i, component in enumerate(components):
                            if complement[i] == "0":
                                continue
                            if product == "map02":
                                # Filter-and-bin map
                                fname_in = os.path.join(
                                    rootdir,
                                    f"{component}_sim/outputs_rk/coadd/{alt_telescope}",
                                    f"coadd_{alt_telescope}_{alt_band}_"
                                    f"{i_time_split:03}of{n_time_split:03}_map.fits",
                                )
                            else:
                                msg = f"Don't know how to assemble signal product: {product}"
                                raise RuntimeError(msg)
                            print(f"Loading {fname_in}", flush=True)
                            m = hp.read_map(fname_in, None)
                            if total is None:
                                total = m
                                component_names = f"{component}"
                            else:
                                total += m
                                component_names += f", {component}"
                        column_names = ["TEMPERATURE", "Q_POLARIZATION", "U_POLARIZATION"]
                        column_units = "K_CMB"
                        signal_header = [
                            ("HIERARCH complement", complement, "Component bit mask"),
                            ("HIERARCH components", component_names, "Component names"),
                        ]
                        print(f"Writing {fname_out}", flush=True)
                        hp.write_map(
                            fname_out,
                            total,
                            nest=False,
                            dtype=np.float32,
                            coord="C",
                            column_names=column_names,
                            column_units=column_units,
                            extra_header=supporting_header + signal_header,
                        )
                        sys.exit()
