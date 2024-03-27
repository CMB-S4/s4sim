import os
import sys

import astropy.units as u
import healpy as hp
import matplotlib.pyplot as plt
#from mpi4py import MPI
import numpy as np
# Numba is not available for Python 3.11 ...
from numba import jit

import toast
from toast.ops import LoadHDF5

#comm = MPI.COMM_WORLD
#rank = comm.rank
#ntask = comm.size
comm = None
rank = 0
ntask = 1

"""
chile_schedule_lat.pruned.txt:
...
 2027-07-02 00:04:59  2027-07-02 01:04:59    61588.003472   61588.045139   180.00 SETTING_SCAN_40                       210.00   330.00    40.00 S       -28.47   284.73   -41.93   281.03   -58.90   283.71   -71.53   285.31  0.07   156   0
 2027-07-02 01:04:59  2027-07-02 02:04:59    61588.045139   61588.086806   180.00 SETTING_SCAN_40                       210.00   330.00    40.00 S       -41.93   281.03   -55.56   277.73   -71.53   285.31   -83.56   305.51  0.06   156   1
 2027-07-02 02:04:59  2027-07-02 03:04:59    61588.086806   61588.128472   180.00 SETTING_SCAN_40                       210.00   330.00    40.00 S       -55.56   277.73   -69.28   274.72   -83.56   305.51   -81.65    62.11  0.06   156   2
 2027-07-02 03:04:59  2027-07-02 04:04:59    61588.128472   61588.170139   180.00 SETTING_SCAN_40                       210.00   330.00    40.00 S       -69.28   274.72   -83.07   272.26   -81.65    62.11   -69.38    74.87  0.06   156   3
 2027-07-02 04:04:59  2027-07-02 05:04:59    61588.170139   61588.211806   180.00 SETTING_SCAN_40                       210.00   330.00    40.00 S       -83.07   272.26   -83.13    87.77   -69.38    74.87   -56.77    75.66  0.06   156   4
 2027-07-02 05:04:59  2027-07-02 06:04:59    61588.211806   61588.253472   180.00 SETTING_SCAN_40                       210.00   330.00    40.00 S       -83.13    87.77   -69.34    85.32   -56.77    75.66   -44.17    74.07  0.05   156   5
 2027-07-02 06:04:59  2027-07-02 07:04:59    61588.253472   61588.295139   180.00 SETTING_SCAN_40                       210.00   330.00    40.00 S       -69.34    85.32   -55.61    82.31   -44.17    74.07   -31.69    71.37  0.05   156   6
 2027-07-02 07:04:59  2027-07-02 08:04:59    61588.295139   61588.336806   180.00 SETTING_SCAN_40                       210.00   330.00    40.00 S       -55.61    82.31   -41.99    79.01   -31.69    71.37   -19.40    67.83  0.05   156   7
 2027-07-02 08:04:59  2027-07-02 09:04:59    61588.336806   61588.378472   180.00 SETTING_SCAN_40                       210.00   330.00    40.00 S       -41.99    79.01   -28.52    75.31   -19.40    67.83    -7.40    63.38  0.05   156   8
 2027-07-02 09:04:59  2027-07-02 10:04:59    61588.378472   61588.420139   180.00 SETTING_SCAN_40                       210.00   330.00    40.00 S       -28.52    75.31   -15.30    71.05    -7.40    63.38     4.19    57.80  0.05   156   9
 2027-07-02 10:04:59  2027-07-02 11:04:59    61588.420139   61588.461806   180.00 SETTING_SCAN_40                       210.00   330.00    40.00 S       -15.30    71.05    -1.79    65.97     4.19    57.80    14.87    50.72  0.04   156  10
 2027-07-02 11:04:59  2027-07-02 12:04:59    61588.461806   61588.503472   180.00 SETTING_SCAN_40                       210.00   330.00    40.00 S        -1.79    65.97     9.89    59.69    14.87    50.72    24.45    41.62  0.04   156  11
 2027-07-02 12:04:59  2027-07-02 13:04:59    61588.503472   61588.545139   180.00 SETTING_SCAN_40                       210.00   330.00    40.00 S         9.89    59.69    21.29    51.63    24.45    41.62    32.25    29.96  0.04   156  12
 2027-07-02 13:04:59  2027-07-02 14:04:59    61588.545139   61588.586806   180.00 SETTING_SCAN_40                       210.00   330.00    40.00 S        21.29    51.63    31.30    41.02    32.25    29.96    37.42    15.54  0.04   156  13
 2027-07-02 14:04:59  2027-07-02 15:04:59    61588.586806   61588.628472   180.00 SETTING_SCAN_40                       210.00   330.00    40.00 S        31.30    41.02    39.07    27.04    37.42    15.54    39.13   359.23  0.04   156  14
 2027-07-02 15:04:59  2027-07-02 15:54:59    61588.628472   61588.663194   180.00 SETTING_SCAN_40                       210.00   330.00    40.00 S        39.07    27.04    43.00    12.69    39.13   359.23    37.64   345.63  0.03   156  15
 2027-07-02 17:10:00  2027-07-02 18:10:00    61588.715278   61588.756944   180.00 RISING_SCAN_40                         30.00   150.00    40.00 R        43.27   348.99    38.59   331.76    30.98   327.83    22.80   316.69  0.03   157   0
 2027-07-02 18:10:00  2027-07-02 19:09:59    61588.756944   61588.798611   180.00 RISING_SCAN_40                         30.00   150.00    40.00 R        38.59   331.76    30.62   318.05    22.80   316.69    12.98   308.03  0.03   157   1
 2027-07-02 19:09:59  2027-07-02 20:09:59    61588.798611   61588.840278   180.00 RISING_SCAN_40                         30.00   150.00    40.00 R        30.62   318.05    20.48   307.66    12.98   308.03     2.20   301.28  0.03   157   2
 2027-07-02 20:09:59  2027-07-02 21:09:59    61588.840278   61588.881944   180.00 RISING_SCAN_40                         30.00   150.00    40.00 R        20.48   307.66     9.01   299.76     2.20   301.28    -9.55   295.96  0.03   157   3
 2027-07-02 21:09:59  2027-07-02 22:09:59    61588.881944   61588.923611   180.00 RISING_SCAN_40                         30.00   150.00    40.00 R         9.01   299.76    -2.58   293.58    -9.55   295.96   -21.57   291.73  0.02   157   4
 2027-07-02 22:09:59  2027-07-02 23:09:59    61588.923611   61588.965278   180.00 RISING_SCAN_40                         30.00   150.00    40.00 R        -2.58   293.58   -16.28   288.56   -21.57   291.73   -33.88   288.40  0.02   157   5
 2027-07-02 23:09:59  2027-07-03 00:09:59    61588.965278   61589.006944   180.00 RISING_SCAN_40                         30.00   150.00    40.00 R       -16.28   288.56   -29.52   284.33   -33.88   288.40   -46.36   285.98  0.02   157   6
...
"""

inroot = "/global/cfs/cdirs/cmbs4/dc/dc1/staging/noise_sim/outputs_rk/LAT0_CHLAT/f150/"
indir = "SETTING_SCAN_40-156-12"
inpath = os.path.join(inroot, indir)
data = toast.Data()

print(f"Loading data from {inpath}", flush=True)
loader = LoadHDF5(
    volume=inpath,
    detdata=["noexistent entry"],
    shared=["times", "flags", "azimuth"],
    pattern="^(RISING|SETTING).*h5",
    #files=[
    #    "/global/cfs/cdirs/cmbs4/dc/dc1/staging/noise_sim/outputs_rk/LAT0_CHLAT/f150/SETTING_SCAN_40-156-0/SETTING_SCAN_40-156-0_LT22_2776490218.h5",
    #    "/global/cfs/cdirs/cmbs4/dc/dc1/staging/noise_sim/outputs_rk/LAT0_CHLAT/f150/SETTING_SCAN_40-156-0/SETTING_SCAN_40-156-0_LT73_130145039.h5",
    #]
)
loader.exec(data)

obs = data.obs[0]
t = obs.shared["times"].data.copy()
t -= t[0]
az = obs.shared["azimuth"].data.copy()

ind = slice(0, 200000)
t = t[ind]
az = np.degrees(az[ind])

nrow, ncol = 2, 2
fig = plt.figure(figsize=[ncol * 8, nrow * 6])
ax = fig.add_subplot(nrow, ncol, 1)
ax.plot(t, az, '.')
ax.set_xlabel("Time [s]")
ax.set_ylabel("Azimuth [deg]")
ax.set_title("Mount azimuth vs. time")

ax = fig.add_subplot(nrow, ncol, 2)
az_rate = np.diff(az) / np.diff(t)
ax.plot(t[1:], az_rate, '.')
ax.set_xlabel("Time [s]")
ax.set_ylabel("Azimuthal rate [deg/s]")
ax.set_title("Mount azimuthal rate vs. time")

ax = fig.add_subplot(nrow, ncol, 3)
az_accel = np.diff(az_rate) / np.diff(t[1:])
ax.plot(t[2:], az_accel, '.')
ax.set_xlabel("Time [s]")
ax.set_ylabel("Azimuthal acceleration [deg/s^2]")
ax.set_title("Mount azimuthal acceleration vs. time")

ax = fig.add_subplot(nrow, ncol, 4)
ax.plot(az[1:], az_rate, '.')
ax.set_xlabel("Azimuth [deg]")
ax.set_ylabel("Azimuthal rate [deg/s]")
ax.set_title("Mount azimuthal rate vs. azimuth")

fig.tight_layout()
fig.savefig("azimuth.png")
