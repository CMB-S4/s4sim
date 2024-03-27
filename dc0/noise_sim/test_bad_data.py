import glob
import os
import pickle
import sys

from toast.ops.polyfilter.kernels import filter_polynomial

fname = "filter_data_POLE3-104-2_LT237_6_0_2.pck"
with open(fname, "rb") as f:
    order, flags, signals, starts, stops = pickle.load(f)

print(f"Filtering {fname} ... ", end="", flush=True)
filter_polynomial(order, flags, signals, starts, stops)
print(f" done!", flush=True)
