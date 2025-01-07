import os
import sys
import warnings

from erfa import ErfaWarning
import astropy.units as u
from astropy.time import Time
import numpy as np
import toml


# astropy.time.Time throws ErfaWarnings for future dates
warnings.filterwarnings("ignore", category=ErfaWarning)

fname = "catalog.txt"

catalog = {}

# Example static source

catalog["static_source"] = {
    "ra_deg" : 30,
    "dec_deg" : -30,
    "freqs_ghz" : [1., 1000.],
    "flux_density_Jy" : [10., 1.],
    "pol_frac" : 0.1,
    "pol_angle_deg" : 0,
}

# Example variable source (timestamps cover the entire simulation span)

catalog["variable_source"] = {
    "ra_deg" : 30,
    "dec_deg" : -25,
    "freqs_ghz" : [1., 1000.],
    "flux_density_Jy" : [
        [10., 1.],
        [30., 10.],
        [10., 1.],
    ],
    "times_mjd" : Time(['2030-01-01', '2030-07-01', '2031-01-01']).mjd,
    "pol_frac" : [0.05, 0.15, 0.05],
    "pol_angle_deg" : [45, 45, 45],
}

# Example transient source, the operator will not extrapolate beyond the provided time stamps

catalog["transient_source"] = {
    "ra_deg" : 30,
    "dec_deg" : -20,
    "freqs_ghz" : [1., 1000.],
    "flux_density_Jy" : [
        [10., 1.],
        [30., 10.],
        [10., 1.],
    ],
    "times_mjd" : Time(
        ['2030-06-14T00:00:00', '2030-06-15T00:00:00', '2030-06-16T00:00:00']
    ).mjd,
}

with open(fname, "w") as f:
    f.write(toml.dumps(catalog))

# Test loading the catalog

with open(fname, "r") as f:
    catalog2 = toml.loads(f.read())
