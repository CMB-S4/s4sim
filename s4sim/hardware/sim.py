# Copyright (c) 2020-2023 CMB-S4 Collaboration.
# Full license can be found in the top level "LICENSE" file.
"""Focalplane simulation tools.
"""

import re
from collections import OrderedDict

import astropy.units as u
import numpy as np
from toast.instrument_coords import xieta_to_quat, quat_to_xieta
from toast.instrument_sim import (
    hex_nring,
    hex_xieta_row_col,
    hex_layout,
    rhomb_dim,
    rhomb_xieta_row_col,
    rhombus_layout,
    rhomb_gamma_angles_qu,
)
import toast.qarray as qa


def ang_to_quat(offsets):
    """Convert cartesian angle offsets and rotation into quaternions.

    Each offset contains two angles specifying the distance from the Z axis
    in orthogonal directions (called "X" and "Y").  The third angle is the
    rotation about the Z axis.  A quaternion is computed that first rotates
    about the Z axis and then rotates this axis to the specified X/Y angle
    location.

    Args:
        offsets (list of arrays):  Each item of the list has 3 elements for
            the X / Y angle offsets in radians and the rotation in radians
            about the Z axis.

    Returns:
        (list): List of quaternions, one for each item in the input list.

    """
    out = list()

    zaxis = np.array([0, 0, 1], dtype=np.float64)

    for off in offsets:
        angrot = qa.rotation(zaxis, off[2])
        wx = np.sin(off[0])
        wy = np.sin(off[1])
        wz = np.sqrt(1.0 - (wx * wx + wy * wy))
        wdir = np.array([wx, wy, wz])
        posrot = qa.from_vectors(zaxis, wdir)
        out.append(qa.mult(posrot, angrot))

    return out


def sim_detectors_toast(hw, tele, tube_slots=None):
    """Update hardware model with simulated detector positions.

    Given a Hardware model, generate all detector properties for the specified
    telescope and optionally a subset of optics tube slots (for the LAT).  The
    detector dictionary of the hardware model is updated in place.
    This function requires the toast subpackage (and hence toast) to be
    importable.
    Args:
        hw (Hardware): The hardware object to update.
        tele (str): The telescope name.
        tube_slots (list, optional): The optional list of tube slots to include.

        Returns:
        None
    """
    sim_telescope_detectors(
        hw, tele, tube_slots=tube_slots,
    )


def sim_detectors_physical_optics(hw, tele, tube_slots=None):
    """Update hardware model with simulated detector positions.
    Given a Hardware model, generate all detector properties for the specified
    telescope and optionally a subset of optics tube slots (for the LAT).  The
    detector dictionary of the hardware model is updated in place.
    This function uses information from physical optics simulations to estimate
    the location of detectors.
    Args:
        hw (Hardware): The hardware object to update.
        tele (str): The telescope name.
        tube_slots (list, optional): The optional list of tube slots to include.

    Returns:
        None
    """
    raise NotImplementedError("Not yet implemented")


def rhombus_hex_layout(
    rhombus_npos, rhombus_width, gap, rhombus_rotate=None, killpix=None
):
    """
    Construct a hexagon from 3 rhombi.

    Args:
        rhombus_npos (int): The number of positions in one rhombus.
        rhombus_width (float): The angle (in degrees) subtended by the
            width of one rhombus along the X axis.
        gap (float): The gap between the edges of the rhombi, in degrees.
        rhombus_rotate (array, optional): An additional angle rotation of
            each position on each rhombus before the rhombus is rotated
            into place.
        killpix (list, optional): Pixel indices to remove for mechanical
            reasons.

    Returns:
        (dict): Keys are the hexagon position and values are quaternions.

    """

    # rhombus dim
    dim = rhomb_dim(rhombus_npos)

    # width in radians
    width_rad = rhombus_width.to_value(u.radian)

    # First layout one rhombus
    rquat = rhombus_layout(
        rhombus_npos,
        rhombus_width,
        "",
        "",
        rhombus_rotate,
    )

    # angular separation of rhombi
    gap_rad = gap.to_value(u.radian)

    # Quaternion offsets of the 3 rhombi
    qcenters = [
        xieta_to_quat(
            0.25 * np.sqrt(3.0) * width_rad + 0.5 * gap_rad,
            -0.25 * width_rad - gap_rad / (2 * np.sqrt(3.0)),
            np.pi / 6,
        ),
        xieta_to_quat(
            0.0,
            0.5 * width_rad + gap_rad / np.sqrt(3.0),
            -0.5 * np.pi,
        ),
        xieta_to_quat(
            -0.25 * np.sqrt(3.0) * width_rad - 0.5 * gap_rad,
            -0.25 * width_rad - gap_rad / (2 * np.sqrt(3.0)),
            5 * np.pi / 6,
        ),
    ]

    nkill = len(killpix)
    result = np.zeros((3 * rhombus_npos - nkill, 4), dtype=np.float64)

    off = 0
    px = 0
    ndigit = int(np.log10(rhombus_npos)) + 1
    nameformat = "{{:0{}}}".format(ndigit)
    for qc in qcenters:
        for p in range(rhombus_npos):
            if px not in killpix:
                name = nameformat.format(p)
                result[off] = qa.mult(qc, rquat[name]["quat"])
                off += 1
            px += 1

    return result


def sector(npos, pos):
    """Return the sector of a given position.
        
        For a hexagonal layout, indexed in a "spiral" scheme (see hex_layout),
        this function returns the sector of a single position.
        The row is zero along the main vertex-vertex axis, and is positive
        or negative above / below this line of positions.
        
        Args:
        npos (int): The number of positions.
        pos (int): position number
        
        Returns:
        (integer): sector of the position
        
        """
    test = npos - 1
    nrings = 1
    while (test - 6 * nrings) >= 0:
        test -= 6 * nrings
        nrings += 1
    if pos == 0:
        sector = 3
    else:
        test = pos - 1
        ring = 1
        while (test - 6 * ring) >= 0:
            test -= 6 * ring
            ring += 1
        if test == 0:
            sector = 3
        else:
            sector = int(test / ring)
    return sector


def sector2(npos, pos):
    """Return the sector of a given position.
        
        For a hexagonal layout, indexed in a "spiral" scheme (see hex_layout),
        this function returns the sector of a single position.
        The row is zero along the main vertex-vertex axis, and is positive
        or negative above / below this line of positions.
        
        Args:
        npos (int): The number of positions.
        pos (int): position number
        
        Returns:
        (integer): sector of the position
        
        """
    test = npos - 1
    nrings = 1
    while (test - 6 * nrings) >= 0:
        test -= 6 * nrings
        nrings += 1
    if pos == 0:
        sector = 0
    else:
        test = pos - 1
        ring = 1
        while (test - 6 * ring) >= 0:
            test -= 6 * ring
            ring += 1
        if test == 0:
            sector = 0
        elif test == ring:
            sector = 0
        else:
            sector = int(test / ring)
    return sector


def triangle(npos, width, rotate=None):
    """Compute positions in an equilateral triangle layout.
        
        Args:
        npos (int): The number of positions packed onto wafer=3
        width (float): distance between tubes in degrees
        rotate (array, optional): Optional array of rotation angles in degrees
        to apply to each position.
        
        Returns:
        (array): Array of quaternions for the positions.
        
        """
    zaxis = np.array([0, 0, 1], dtype=np.float64)
    sixty = np.pi / 3.0
    thirty = np.pi / 6.0
    rtthree = np.sqrt(3.0)
    rtthreebytwo = 0.5 * rtthree

    tubedist = width * np.pi / 180.0
    result = np.zeros((npos, 4), dtype=np.float64)
    posangarr = np.array([sixty * 3.0 + thirty, -thirty, thirty * 3.0])
    for pos in range(npos):
        posang = posangarr[pos]
        posdist = tubedist / rtthree

        posx = np.sin(posdist) * np.cos(posang)
        posy = np.sin(posdist) * np.sin(posang)
        posz = np.cos(posdist)
        posdir = np.array([posx, posy, posz], dtype=np.float64)
        norm = np.sqrt(np.dot(posdir, posdir))
        posdir /= norm
        posrot = qa.from_vectors(zaxis, posdir)

        if rotate is None:
            result[pos] = posrot
        else:
            prerot = qa.rotation(zaxis, rotate[pos] * np.pi / 180.0)
            result[pos] = qa.mult(posrot, prerot)

    return result


def kill_pixels(layout, killpix):
    if killpix is None or len(killpix) == 0:
        return
    # Kill pixels
    for pix in killpix:
        name = f"{pix:03}"
        del layout[name]
    # Assign new pixel numbers
    for offset, old_name in enumerate(sorted(layout.keys())):
        new_name = f"{offset:03}"
        if old_name == new_name:
            continue
        layout[new_name] = layout[old_name]
        del layout[old_name]
    return

def sim_wafer_detectors(
    hw,
    wafer,
    platescale,
    fwhm,
    band=None,
    partial_type=None,
    no_gap=None,
    center=np.array([0, 0, 0, 1], dtype=np.float64),
):
    """Generate detector properties for a wafer.

    Given a Hardware configuration, generate all detector properties for
    the specified wafer and optionally only the specified band.

    Args:
        hw (Hardware): The hardware properties.
        wafer (str): The wafer name.
        platescale (float): The plate scale in degrees / mm.
        fwhm (dict): Dictionary of nominal FWHM values in arcminutes for
            each band.
        band (str, optional): Optionally only use this band.
        center (array, optional): The quaternion offset of the center.

    Returns:
        (OrderedDict): The properties of all selected detectors.

    """
    # The properties of this wafer
    wprops = hw.data["wafers"][wafer]
    # The readout card and its properties
    card = wprops["card"]
    cardprops = hw.data["cards"][card]
    # The bands
    bands = wprops["bands"]
    if band is not None:
        if band in bands:
            bands = [band]
        else:
            raise RuntimeError("band '{}' not valid for wafer '{}'".format(band, wafer))

    # Lay out the pixel locations depending on the wafer type.  Also
    # compute the polarization orientation rotation, as well as the A/B
    # handedness for the Sinuous detectors.

    npix = wprops["npixel"]
    pixsep = platescale * wprops["pixsize"]
    layout_A = None
    layout_B = None
    handed = None
    if wprops["packing"] == "RP":
        # Feedhorn (NIST style)
        # Make gap zero for partial_arrays
        if no_gap is None:
            gap = platescale * wprops["rhombusgap"]
        else:
            gap = 0.0

        nrhombus = npix // 3
        # This dim is also the number of pixels along the short axis.
        dim = rhomb_dim(nrhombus)
        # This is the center-center distance along the short axis
        width = (dim - 1) * pixsep
        # The orientation within each rhombus alternates between zero and 45
        # degrees.  However there is an offset.  We choose this arbitrarily
        # for the nominal rhombus position, and then the rotation of the
        # other 2 rhombi will naturally modulate this.
        pol_A = np.zeros(nrhombus, dtype=np.float64)
        pol_B = np.zeros(nrhombus, dtype=np.float64)
        poloff = 22.5
        for p in range(nrhombus):
            # get the row / col of the pixel
            row, col = rhomb_xieta_row_col(nrhombus, p)
            if np.mod(row, 2) == 0:
                pol_A[p] = 0.0 + poloff
            else:
                pol_A[p] = 45.0 + poloff
            pol_B[p] = 90.0 + pol_A[p]
        # kill pixels in partial arrays
        if partial_type is None:
            kill = []
        elif partial_type == "rhombus":
            kill = np.linspace(dim ** 2, 3 * dim ** 2 - 1, 2 * dim ** 2, dtype=np.int)
        elif partial_type == "half":
            kill1 = np.linspace(
                0, ((dim ** 2 - dim) // 2 - 1), (dim ** 2 - dim) // 2, dtype=np.int
            )
            kill2 = np.linspace(dim ** 2, 2 * dim ** 2 - 1, dim ** 2, dtype=np.int)
            kill = np.append(kill1, kill2)
        else:
            kill = []
        # We are going to remove 2 pixels for mechanical reasons
        # kf = dim * (dim - 1) // 2
        # kill = [kf, kf + dim - 2]
        layout_A = rhombus_hex_layout(
            nrhombus,
            width * u.degree,
            gap * u.degree,
            rhombus_rotate=pol_A * u.degree,
            killpix=kill,
        )
        layout_B = rhombus_hex_layout(
            nrhombus,
            width * u.degree,
            gap * u.degree,
            rhombus_rotate=pol_B * u.degree,
            killpix=kill,
        )
    elif wprops["packing"] == "HP":
        # Hex close-packed
        # This is the center-center distance along the vertex-vertex axis
        width = (2 * (hex_nring(npix) - 1)) * pixsep
        # The sinuous handedness is chosen so that A/B pairs of pixels have the
        # same nominal orientation but trail each other along the
        # vertex-vertex axis of the hexagon.  The polarization orientation
        # changes every other column
        pol_A = np.zeros(npix, dtype=np.float64)
        pol_B = np.zeros(npix, dtype=np.float64)
        for p in range(npix):
            row, col = hex_xieta_row_col(npix, p)
            if np.mod(col, 4) < 2:
                pol_A[p] = 0.0
            else:
                pol_A[p] = 45.0
            pol_B[p] = 90.0 + pol_A[p]
        if partial_type is None:
            kill = []
        elif partial_type == "half":
            kill = []
            for ii in range(npix):
                if sector(npix, ii) < 3:
                    kill.append(ii)
        elif partial_type == "rhombus":
            kill = []
            for ii in range(npix):
                if (
                    (sector2(npix, ii) == 1)
                    or (sector2(npix, ii) == 2)
                    or (sector2(npix, ii) == 3)
                    or (sector2(npix, ii) == 4)
                ):
                    kill.append(ii)
        else:
            kill = []
        layout_A = hex_layout(
            npix,
            width * u.degree,
            "",
            "",
            pol_A * u.degree,
        )
        kill_pixels(layout_A, kill)
        layout_B = hex_layout(
            npix,
            width * u.degree,
            "",
            "",
            pol_B * u.degree,
        )
        kill_pixels(layout_B, kill)
    else:
        msg = f"Unknown wafer packing '{wprops['packing']}'"
        raise RuntimeError(msg)

    # Now we go through each pixel and create the orthogonal detectors for
    # each band.
    dets = OrderedDict()

    chan_per_coax = cardprops["nchannel"] // cardprops["ncoax"]
    chan_per_bias = cardprops["nchannel"] // cardprops["nbias"]

    doff = 0
    p = 0
    idoff = int(wafer) * 10000
    for px in range(npix):
        if px in kill:
            continue
        pstr = f"{p:03}"
        for b in bands:
            for pl, layout in zip(["A", "B"], [layout_A, layout_B]):
                dprops = OrderedDict()
                dprops["wafer"] = wafer
                dprops["ID"] = idoff + doff
                dprops["pixel"] = pstr
                dprops["band"] = b
                dprops["fwhm"] = fwhm[b]
                dprops["pol"] = pl
                if handed is not None:
                    dprops["handed"] = handed[p]
                # Made-up assignment to readout channels
                dprops["card"] = card
                dprops["channel"] = doff
                dprops["coax"] = doff // chan_per_coax
                dprops["bias"] = doff // chan_per_bias
                # Layout quaternion offset is from the origin.  Now we apply
                # the rotation of the wafer center.
                try:
                    dprops["quat"] = qa.mult(center, layout[p]).flatten()
                except KeyError:
                    dprops["quat"] = qa.mult(center, layout[f"{p:03}"]["quat"]).flatten()
                dname = f"{wafer}_{pstr}_{b}_{pl}"
                dets[dname] = dprops
                doff += 1
        p += 1

    return dets


def sim_telescope_detectors(hw, tele, tubes=None):
    """Generate detector properties for a telescope.

    Given a Hardware model, generate all detector properties for the specified
    telescope and optionally a subset of optics tubes (for the LAT).

    Args:
        hw (Hardware): The hardware object to use.
        tele (str): The telescope name.
        tubes (list, optional): The optional list of tubes to include.

    Returns:
        (OrderedDict): The properties of all selected detectors.

    """
    zaxis = np.array([0, 0, 1], dtype=np.float64)
    thirty = np.pi / 6.0
    # The properties of this telescope
    teleprops = hw.data["telescopes"][tele]
    tele_platescale = teleprops["platescale"]
    fwhm = teleprops["fwhm"]

    # The tubes
    alltubes = teleprops["tubes"]
    ntube = len(alltubes)
    if tubes is None:
        tubes = alltubes
    else:
        for t in tubes:
            if t not in alltubes:
                msg = f"Invalid tube '{t}' for telescope '{tele}'"
                raise RuntimeError(msg)

    alldets = OrderedDict()
    if ntube == 3:
        # This is a SAT.  We have three tubes.
        tubespace = teleprops["tubespace"]
        # tuberot = 0.0 * np.ones(7, dtype=np.float64)
        # tcenters = hex_layout(
        #     7, 2 * (tubespace * tele_platescale), rotate=tuberot
        # )
        # tuberot = 90.0 * np.ones(3, dtype=np.float64)
        # tcenters = triangle(3, (tubespace * tele_platescale), rotate=tuberot)

        for tindx, tube in enumerate(tubes):
            tubeprops = hw.data["tubes"][tube]
            waferspace = tubeprops["waferspace"]
            platescale = tubeprops["platescale"]
            location = str(tubeprops["toast_hex_pos"])
            type = tubeprops["type"]
            if type == "SAT_HF":
                tcenters = hex_layout(
                    7,
                    2 * tubespace * tele_platescale * u.degree,
                    "",
                    "",
                    90 * np.ones(7, dtype=float) * u.degree,
                )
                nwafer = len(tubeprops["wafers"])
                wcenters = hex_layout(
                    nwafer,
                    2 * waferspace * platescale * u.degree,
                    "",
                    "",
                    np.zeros(nwafer, dtype=float) * u.degree,
                )
                centers = list()
                for p, q in wcenters.items():
                    centers.append(qa.mult(tcenters[location]["quat"], q["quat"]))

                for windx, wafer in enumerate(tubeprops["wafers"]):
                    if windx == 4:
                        partial_type = "half"
                    elif windx == 5:
                        partial_type = "half"
                    elif windx == 7:
                        partial_type = "half"
                    elif windx == 8:
                        partial_type = "half"
                    elif windx == 10:
                        partial_type = "half"
                    elif windx == 11:
                        partial_type = "half"
                    else:
                        partial_type = None
                    dets = sim_wafer_detectors(
                        hw,
                        wafer,
                        platescale,
                        fwhm,
                        center=centers[windx],
                        partial_type=partial_type,
                    )
                    alldets.update(dets)
            else:
                tcenters = hex_layout(
                    7,
                    2 * tubespace * tele_platescale * u.degree,
                    "",
                    "",
                    np.zeros(7, dtype=float) * u.degree,
                )
                nwafer = len(tubeprops["wafers"])
                wcenters = hex_layout(
                    nwafer,
                    2 * waferspace * platescale * u.degree,
                    "",
                    "",
                    np.zeros(nwafer, dtype=float) * u.degree,
                )
                centers = list()
                for p, q in wcenters.items():
                    centers.append(qa.mult(tcenters[location]["quat"], q["quat"]))

                for windx, wafer in enumerate(tubeprops["wafers"]):
                    partial_type = None
                    dets = sim_wafer_detectors(
                        hw,
                        wafer,
                        platescale,
                        fwhm,
                        center=centers[windx],
                        partial_type=partial_type,
                    )
                    alldets.update(dets)
    else:
        # This is the LAT.  Compute the tube centers.
        # Rotate each tube by 90 degrees, so that it is pointed "down".
        tubespace = teleprops["tubespace"]
        tcenters = hex_layout(
            91,
            10 * tubespace * tele_platescale * u.degree,
            "",
            "",
            90 * np.ones(91, dtype=float) * u.degree,
        )

        for tindx, tube in enumerate(tubes):
            tubeprops = hw.data["tubes"][tube]
            waferspace = tubeprops["waferspace"]
            platescale = tubeprops["platescale"]
            location = f"{tubeprops['toast_hex_pos']:02}"

            # get centers and rotations for arrays
            wcenters = [np.array([0.0, 0.0, 0.0])]
            qwcenters = ang_to_quat(wcenters)
            centers = list()
            for qwc in qwcenters:
                centers.append(qa.mult(tcenters[location]["quat"], qwc))

            windx = 0
            for wafer in tubeprops["wafers"]:
                # For first three wafers, use whole wafers,
                # then construct partial wafers
                dets = sim_wafer_detectors(
                    hw,
                    wafer,
                    platescale,
                    fwhm,
                    center=centers[windx],
                    partial_type=None,
                    no_gap=None,
                )
                alldets.update(dets)
                #windx += 1
    return alldets
