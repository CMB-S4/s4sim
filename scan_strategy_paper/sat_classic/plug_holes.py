# This file contains one utility routine lifted
# from toast-npipe because toast-npipe requires
# a functional toast-2 install to work.

import os
import sys

import healpy as hp
import numpy as np


def plug_holes(m, verbose=False, in_place=True, nest=False):
    """
    Use simple downgrading to derive estimates of the missing pixel values
    """
    nbad_start = np.sum(np.isclose(m, hp.UNSEEN))

    if nbad_start == m.size:
        if verbose:
            print("plug_holes: All map pixels are empty. Cannot plug holes", flush=True)
        return

    if nbad_start == 0:
        return

    nside = hp.get_nside(m)
    npix = m.size
    if nest:
        mnest = m.copy()
    else:
        mnest = hp.reorder(m, r2n=True)

    lowres = mnest
    nside_lowres = nside
    bad = np.isclose(mnest, hp.UNSEEN)
    while np.any(bad) and nside_lowres > 1:
        nside_lowres //= 2
        lowres = hp.ud_grade(lowres, nside_lowres, order_in="NESTED")
        hires = hp.ud_grade(lowres, nside, order_in="NESTED")
        bad = np.isclose(mnest, hp.UNSEEN)
        mnest[bad] = hires[bad]

    nbad_end = np.sum(bad)

    if nbad_end != 0:
        mn = np.mean(mnest[np.logical_not(bad)])
        mnest[bad] = mn

    if not in_place:
        m = m.copy()
    if nest:
        m[:] = mnest
    else:
        m[:] = hp.reorder(mnest, n2r=True)

    if verbose and nbad_start != 0:
        print(
            "plug_holes: Filled {} missing pixels ({:.2f}%), lowest "
            "resolution was Nside={}.".format(
                nbad_start, (100.0 * nbad_start) // npix, nside_lowres
            )
        )
    return m
