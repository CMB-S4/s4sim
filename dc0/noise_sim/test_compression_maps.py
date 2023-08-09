import os
import sys

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np

from toast.pixels_io_healpix import read_healpix, write_healpix


fname_test = "RISING_SCAN_40-10-7_LT38_2513396920.h5"
fname_map = "MapMaker_binmap.h5"

refdir = "outputs/LAT0_CHLAT/f090/RISING_SCAN_40-10-7"
ref = read_healpix(os.path.join(refdir, fname_map), verbose=False)
refsize = os.path.getsize(os.path.join(refdir, fname_test))

for precision in range(2, 9):
    testdir = f"testdata_precision{precision}"
    ver = read_healpix(os.path.join(testdir, fname_map), verbose=False)
    versize = os.path.getsize(os.path.join(testdir, fname_test))

    print("")
    print(f"Precision = {precision}")
    print(f"Compression ratio = {refsize / versize:.3f}, file size = {100 * versize / refsize:.3f}%")
    print(f"{refsize / versize:.3f}   {100 * versize / refsize:.3f}%  ", end="")
    #print(f"Relative compression error:", end="")
    for i in range(3):
        good = ref[i] != 0
        rms0 = np.std(ref[i][good])
        rms1 = np.std((ref[i] - ver[i])[good])
        print(f" {rms1 / rms0}", end="")
    print("", flush=True)
    
    del ver
