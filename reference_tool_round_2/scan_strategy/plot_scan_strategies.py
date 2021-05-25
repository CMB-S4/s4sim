import os
import sys

import matplotlib.pyplot as plt
import healpy as hp

for site in "pole", "chile":
    for tele in "sat", "lat":
        fname = "{}_{}/out/00000000/toast_telescope_all_time_all_hmap.fits".format(
            site, tele
        )
        hmap = hp.read_map(fname, verbose=False)
        hmap[hmap == 0] = hp.UNSEEN
        hp.mollview(hmap, title="{} {}".format(site, tele), xsize=800, unit="hits")
        hp.graticule(22.5, verbose=False)
        plt.savefig("hits_{}_{}.png".format(site, tele))
