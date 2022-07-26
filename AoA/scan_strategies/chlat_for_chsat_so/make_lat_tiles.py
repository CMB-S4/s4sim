import numpy as np
import os
import sys
import matplotlib.pyplot as plt

import healpy as hp

plt.style.use('classic')

#dec_min, dec_max, dec_step, dec_overlap = -75, 75, 30, 15
dec_min, dec_max, dec_step, dec_overlap = -80, 80, 20, 10
ra_step, ra_overlap = 10, 5
# Use this power to adjust contrast in the priorities
# map_power = 0

nside = 512
npix = 12 * nside ** 2
lon, lat = hp.pix2ang(nside, np.arange(npix), lonlat=True)

# lon_min, lat_max, lon_max, lat_min, priority

lat_patches = {
    #'NAME' : [lonmin, latmax, lonmax, latmin, priority]
    #'south': [-50, -30, 90, -50],
    'south1': [-60, -20, 100, -40, 0.001],  # The south patch is the primary patch
    'south2': [-70, -30, 110, -50, 0.001],  # The south patch is the primary patch
    'south3': [-60, -40, 100, -60, 0.001],  # The south patch is the primary patch
    #'north': [120,  0, 195, -20]
    'north1': [120, 10, 195, -10, 1000],  # The Northern patch is secondary
    'north2': [110,  0, 205, -20, 1000],  # The Northern patch is secondary
    'north3': [120,-10, 195, -30, 1000],  # The Northern patch is secondary
}

def unwind(alpha, beta):
    # Find the branch in beta that is closest to alpha
    while np.abs(alpha - beta + 360) < np.abs(alpha - beta):
        beta -= 360
    while np.abs(alpha - beta - 360) < np.abs(alpha - beta):
        beta += 360
    return beta

patch_map = np.zeros(npix)
hit_map = np.zeros(npix)

def plot_patches(patches):
    for patch in patches:
        lonmin, latmax, lonmax, latmin, priority = patches[patch]
        lons, lats = [], []
        for lon in np.linspace(lonmin, lonmax, 100):
            lons.append(lon)
            lats.append(latmax)
        for lat in np.linspace(latmax, latmin, 100):
            lons.append(lonmax)
            lats.append(lat)
        for lon in np.linspace(lonmax, lonmin, 100):
            lons.append(lon)
            lats.append(latmin)
        for lat in np.linspace(latmin, latmax, 100):
            lons.append(lonmin)
            lats.append(lat)
        hp.projplot(lons, lats, '-', color='white', lw=2, threshold=1,
                    lonlat=True, coord='C')

with open('patches_lat.txt', 'w') as fout:
    # Add entries for the full patches
    for name in lat_patches:
        lonmin, latmax, lonmax, latmin, priority = lat_patches[name]
        fout.write('--patch\n')
        fout.write('{},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f}\n'.format(
            name, priority, lonmin, latmax, lonmax, latmin))
    dec_start = dec_min
    while dec_start < dec_max:
        dec_stop = min(dec_start + dec_step, dec_max)
        latgood = np.logical_and(lat >= dec_start, lat <= dec_stop)
        # Adjust the RA step by latitude so that the patch areas
        # are more uniform
        if dec_start * dec_stop < 0:
            minlat = 0
        else:
            minlat = min(np.abs(dec_start), np.abs(dec_stop))
        latfac = 1 / np.cos(np.radians(minlat))
        ra_step_lat = ra_step * latfac
        nstep = int(360 / ra_step_lat)
        ra_step_lat = 360 / nstep
        ra_overlap_lat = ra_overlap * ra_step_lat / ra_step
        ra_start = 0
        for ilon in range(nstep):
            # Tile extent in RA
            ra_start = ilon * ra_step_lat
            ra_stop = ra_start + ra_step_lat
            # see if the tile overlaps with the LAT patches
            overlaps = False
            for patch in lat_patches:
                lonmin, latmax, lonmax, latmin, priority = lat_patches[patch]
                if latmax - dec_overlap > dec_start and \
                   latmin + dec_overlap < dec_stop:
                    # Tile overlaps with the patch in DEC.
                    # How about RA?
                    ra1 = unwind(lonmin, ra_start)
                    ra2 = ra_stop + ra1 - ra_start
                     # Accept any overlap
                    if lonmax > ra1 + 1 and lonmin < ra2 - 1:
                        print('RA patch = {:4}:{:4}, '
                              'RA tile = {:4}:{:4}'.format(
                            int(lonmin), int(lonmax), int(ra1), int(ra2)))
                    # Require complete overlap
                    #if lonmin <= ra1 and ra2 <= lonmax:
                        overlaps = True
                        break
            if overlaps:
                name = 'DEC{:+04}..{:+04}_RA{:+08.3f}..{:+08.3f}'.format(
                    dec_start, dec_stop, ra_start, ra_stop)
                good = np.logical_and(lon >= ra_start, lon <= ra_stop) * latgood
                patch_map[good] += 1 / priority
                hit_map[good] += 1
                fout.write('--patch\n')
                fout.write('{},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f}\n'.format(
                    name, priority * 1000, ra_start, dec_stop, ra_stop, dec_start))
            ra_start += ra_step_lat - ra_overlap_lat
        dec_start += dec_step - dec_overlap

good = hit_map != 0
patch_map[good] /= hit_map[good]
patch_map[patch_map == 0] = hp.UNSEEN
hp.mollview(patch_map, title='priority')
hp.graticule(15, verbose=False)
plot_patches(lat_patches)
plt.savefig('priority_lat.png')
