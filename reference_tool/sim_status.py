#!/usr/bin/env python

import os
import sys
from collections import OrderedDict

BLACK = '\033[0;30m'
RED = '\033[0;31m'
GREEN = '\033[0;32m'
ORANGE = '\033[0;33m'
BLUE = '\033[0;34m'
PURPLE = '\033[0;35m'
CYAN = '\033[0;36m'
LGRAY = '\033[0;37m'
DGREY = '\033[1;30m'
LRED = '\033[1;31m'
LGREEN = '\033[1;32m'
YELLOW = '\033[1;33m'
LBLUE = '\033[1;34m'
LPURPLE = '\033[1;35m'
LCYAN = '\033[1;36m'
WHITE = '\033[1;37m'

CLEAR = '\033[0m'
BOLD = '\033[1m'
# FAINT = '\033[2m'
UNDERLINE = '\033[4m'
REVERSE = '\033[7m'
# CROSSED_OUT = '\033[9m'
# FRAMED = '\033[51m'
# ENCIRCLED = '\033[52m'
# OVERLINED = '\033[53m'

flavors = (
    "noise",
    "atmosphere",
    "cmb-unlensed",
    "cmb-lensing",
    "cmb-tensors",
    "foreground",
    #"cmb-scalar",
    #"cmb-lensing",
    #"cmb-tensors",
    #"galactic",
    #"extra-galactic",
)

telescopes = OrderedDict()
telescopes["LAT"] = OrderedDict()

lat = telescopes["LAT"]
lat["LT56"] = ["ULFL1"]  #                  20 GHz
lat["LT17"] = ["LFL1", "LFL2"]  #     27 &  39 GHz
lat["LT5"] = ["MFL1", "MFL2"]  #      93 & 145 GHz
lat["LT0"] = ["HFL1", "HFL2"]  #     225 & 278 GHz
#lat["LT56"] = ["ULFPL1"]  #                 20 GHz Pole
#lat["LT54"] = ["LFPL1", "LFPL2"]  #   27 &  39 GHz Pole
#lat["LT42"] = ["MFPL1", "MFPL2"]  #   93 & 145 GHz Pole
#lat["LT38"] = ["HFPL1", "HFPL2"]  #  225 & 278 GHz Pole

telescopes["SAT"] = OrderedDict()
sat = telescopes["SAT"]
sat["ST16"] = ["LFS1", "LFS2"]  #   30 &  40 GHz - SAT5 - FOV 17.5 deg
sat["ST0"] = ["MFLS1", "MFLS2"]  #  85 & 145.1 GHz - SAT0 - FOV 14.5 deg
sat["ST6"] = ["MFHS1", "MFHS2"]  #  95 & 155.1 GHz - SAT2 - FOV 14.5 deg
sat["ST12"] = ["HFS1", "HFS2"]  #  220 & 270 GHz - SAT4 - FOV 17.5 deg

nmc = 8
rootdir = "out"
print('\nSimulation status:\n')

sites = "chile", "pole"
print("{:12}{:90}{:90}".format("", "Chile", "Pole"))
print("{:12}{:42}{:48}{:42}{:48}".format("", "LAT", "SAT", "LAT", "SAT"))

print("{:12}".format("flavor"), end="")
for site in sites:
    for telescope, tubes in telescopes.items():
        for tube, bands in tubes.items():
            for band in bands:
                #if site == "chile" and telescope == "LAT" and "P" in band and band != "ULFPL1":
                #    continue
                #if site == "chile" and telescope == "LAT" and "P" not in band:
                #    continue
                print("{:6}".format(band), end="")
print()

for flavor in flavors:
    print("{:12}".format(flavor), end="")
    for site in sites:
        for telescope, tubes in telescopes.items():
            for tube, bands in tubes.items():
                for band in bands:
                    #if site == "chile" and telescope == "LAT" and "P" in band and band != "ULFPL1":
                    #    continue
                    #if site == "chile" and telescope == "LAT" and "P" not in band:
                    #    continue
                    for mc in range(100):
                        fname = "out/{:08}/{}_{}_{}_{}_filtered_telescope_all_time_all_bmap.fits".format(mc, site, flavor, telescope, band)
                        if not os.path.isfile(fname):
                            break
                    if mc == 0:
                        print("{:6}".format(""), end="")
                    else:
                        print("{:3}   ".format(mc), end="")
    print()
    """
    mapver = {"nominal" : "nominal", "realistic" : "perturbed"}[ver]
    print('\n{} case:'.format(ver))
    print(' ' * 5, end='')
    for tube in tubes:
        print('{:10}|'.format(tube), end='')
    print("")
    print(' ' * 5, end='')
    for tube in tubes:
        bands = bands_by_tube[tube]
        for band in bands:
            print('{:5}'.format(band), end='')
        print('|', end="")
    print()
    for isub in range(1, nsub + 1):
        if isub == 11:
            print('{}{}'.format('-' * 5, '----------+' * len(tubes)))
        print('{:02}/{:02}'.format(isub, nsub), end="")
        for tube in tubes:
            bands = bands_by_tube[tube]
            for band in bands:
                fname0 = "{}_{}_{:02}_of_{:02}.{}.log".format(tube, band, isub, nsub, ver)
                fname1 = os.path.join(
                    rootdir,
                    "{}_{}_{:02}_of_{:02}.{}_telescope_all_time_all_hmap.fits".format(
                        tube, band, isub, nsub, mapver))
                fname2 = os.path.join(
                    rootdir,
                    "{}_{}_{:02}_of_{:02}.{}_filtered_telescope_all_time_all_bmap.fits".format(
                        tube, band, isub, nsub, mapver))
                res = RED + REVERSE + '     ' + CLEAR
                if os.path.isfile(fname0):
                    res = YELLOW + REVERSE + '  .  ' + CLEAR
                if os.path.isfile(fname1):
                    res = GREEN + REVERSE + '  +  ' + CLEAR
                if os.path.isfile(fname2):
                    res = BLUE + REVERSE + '  X  ' + CLEAR
                print('{:5}'.format(res), end='')
            print('|', end='')
        print()
#print(BOLD + '\nRealizations done: {:.3f}\n'.format(total_completed / 200) + CLEAR)
print('Legend:')
print(RED + REVERSE + '     ' + CLEAR +  ' = not done')
print(YELLOW + REVERSE + '  .  ' + CLEAR + ' = began running')
print(GREEN + REVERSE + '  +  ' + CLEAR + ' = destriped map done')
print(BLUE + REVERSE + '  X  ' + CLEAR + ' = filtered map done')
#print(CYAN + REVERSE + '  O  ' + CLEAR + ' = single detector pol. templates done')

"""
