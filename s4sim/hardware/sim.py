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
    rhombus_hex_layout,
)
import toast.qarray as qa


XAXIS, YAXIS, ZAXIS = np.eye(3)


def sim_detectors_toast(hw, tele, tubes=None):
    """Update hardware model with simulated detector positions.

    Given a Hardware model, generate all detector properties for the specified
    telescope and optionally a subset of optics tube slots (for the LAT).  The
    detector dictionary of the hardware model is updated in place.
    This function requires the toast subpackage (and hence toast) to be
    importable.
    
    Args:
        hw (Hardware): The hardware object to update.
        tele (str): The telescope name.
        tubes (list, optional): The optional list of tube slots to include.

    Returns:
        None
    
    """
    sim_telescope_detectors(
        hw, tele, tubes=tubes,
    )


def sim_detectors_physical_optics(hw, tele, tubes=None):
    """Update hardware model with simulated detector positions.

    Given a Hardware model, generate all detector properties for the specified
    telescope and optionally a subset of optics tube slots (for the LAT).  The
    detector dictionary of the hardware model is updated in place.
    This function uses information from physical optics simulations to estimate
    the location of detectors.
    
    Args:
        hw (Hardware): The hardware object to update.
        tele (str): The telescope name.
        tubes (list, optional): The optional list of tube slots to include.

    Returns:
        None
    
    """
    raise NotImplementedError("Not yet implemented")


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


def sim_wafer_detectors(
    hw,
    wafer,
    platescale,
    fwhm,
    band=None,
    partial_type=None,
    center=np.array([0, 0, 0, 1], dtype=np.float64),
):
    """Generate detector properties for a wafer.

    Given a Hardware configuration, generate all detector properties for
    the specified wafer and optionally only the specified band.

    By convention, polarization angle quantities are stored in degrees and
    fwhm is stored in arcminutes.  These units are stripped to ease serialization
    to TOML and JSON.  The units are restored when constructing focalplane tables
    for TOAST.

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
            msg = f"band '{band}' not valid for wafer '{wafer}'"
            raise RuntimeError(msg)

    # Lay out the pixel locations depending on the wafer type.  Also
    # compute the polarization orientation rotation, as well as the A/B
    # handedness for the Sinuous detectors.

    npix = wprops["npixel"]
    pixsep = platescale * wprops["pixsize"] * u.degree
    handed = None
    if wprops["packing"] == "RP":
        # Feedhorn (NIST style)
        # The "gap" is the **additional** space beyond the normal pixel
        # separation.
        gap = 0.0 * u.degree

        # Orientation of the wafer with respect to the focalplane
        wafer_rot_rad = 0.0

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
                pol_A[p] = 0 + poloff
            else:
                pol_A[p] = 45 + poloff
            pol_B[p] = 90 + pol_A[p]
        pol_A = u.Quantity(pol_A, u.degree)
        pol_B = u.Quantity(pol_B, u.degree)

        layout_A = rhombus_hex_layout(
            nrhombus, width, "", "", gap=gap, pol=pol_A
        )
        layout_B = rhombus_hex_layout(
            nrhombus, width, "", "", gap=gap, pol=pol_B
        )

        # kill pixels in partial arrays
        if partial_type is None:
            kill = wprops["pins"]
        elif partial_type == "rhombus":
            kill = np.linspace(dim**2, 3 * dim**2 - 1, 2 * dim**2, dtype=np.int32)
        elif partial_type == "half":
            kill1 = np.linspace(
                0, ((dim**2 - dim) // 2 - 1), (dim**2 - dim) // 2, dtype=np.int32
            )
            kill2 = np.linspace(dim**2, 2 * dim**2 - 1, dim**2, dtype=np.int32)
            kill = np.append(kill1, kill2)
        else:
            kill = wprops["pins"]
    elif wprops["packing"] == "HP":
        # Hex close-packed
        # This is the center-center distance along the vertex-vertex axis
        width = (2 * (hex_nring(npix) - 1)) * pixsep

        # The orientation of the wafer with respect to the focalplane.
        wafer_rot_rad = 0.0
        hex_cent = qa.from_axisangle(ZAXIS, wafer_rot_rad)

        # The sinuous handedness is chosen so that A/B pairs of pixels have the
        # same nominal orientation but trail each other along the
        # vertex-vertex axis of the hexagon.  The polarization orientation
        # changes every other column
        pol_A = np.zeros(npix, dtype=float)
        pol_B = np.zeros(npix, dtype=float)
        for p in range(npix):
            row, col = hex_xieta_row_col(npix, p)
            if np.mod(col, 4) < 2:
                pol_A[p] = 0.0
            else:
                pol_A[p] = 45.0
            pol_B[p] = 90 + pol_A[p]
        pol_A = u.Quantity(pol_A, u.degree)
        pol_B = u.Quantity(pol_B, u.degree)

        layout_A = hex_layout(npix, width, "", "", pol_A, center=hex_cent)
        layout_B = hex_layout(npix, width, "", "", pol_B, center=hex_cent)
        
        if partial_type is None:
            kill = wprops["pins"]
        elif partial_type == "half":
            kill = []
            for ii in range(npix):
                if sector(npix, ii) < 3:
                    kill.append(ii)
        elif partial_type == "rhombus":
            kill = []
            for ii in range(npix):
                if (sector2(npix, ii) in [1, 2, 3, 4]):
                    kill.append(ii)
        else:
            kill = wprops["pins"]
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
        if npix < 100:
            pxstr = f"{px:02d}"
        else:
            pxstr = f"{px:03d}"
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
                    dprops["handed"] = handed[px]
                # Made-up assignment to readout channels
                dprops["card"] = card
                dprops["channel"] = doff
                dprops["coax"] = doff // chan_per_coax
                dprops["bias"] = doff // chan_per_bias

                # Polarization angle in wafer basis.  This is the gamma angle
                # returned by the layout functions above, less the rotation
                # of the wafer.
                wafer_gamma = layout[pxstr]["gamma"]
                dprops["pol_ang_wafer"] = np.degrees(wafer_gamma - wafer_rot_rad)

                # Polarization angle in focalplane basis.  This is, by
                # definition, equal to the gamma angle computed by the layout
                # functions.  However, we must also account for the extra rotation
                # of the center offset passed to this function.
                if center is not None:
                    dprops["quat"] = qa.mult(center, layout[pxstr]["quat"]).flatten()
                    _, _, temp_gamma = quat_to_xieta(dprops["quat"])
                    dprops["gamma"] = np.degrees(temp_gamma)
                else:
                    dprops["quat"] = layout[pxstr]["quat"].flatten()
                    dprops["gamma"] = np.degrees(wafer_gamma)
                dprops["pol_ang"] = dprops["gamma"]

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

    # Wafer props
    wprops = hw.data["wafers"]

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
        # This is a SAT.  We have three co-incident tubes.
        tubespace = teleprops["tubespace"]

        for tindx, tube in enumerate(tubes):
            tubeprops = hw.data["tubes"][tube]
            waferspace = tubeprops["waferspace"]
            platescale = tubeprops["platescale"]
            location = str(tubeprops["toast_hex_pos"])
            wafer_ang_deg = np.array(tubeprops["wafer_angle"])
            wafer_ang_rad = np.radians(wafer_ang_deg)
            type = tubeprops["type"]

            # We only want 12 wafers but we want to offset the
            # default 19-wafer hexagonal layout by half a wafer
            # and then drop the ones which are the furthest away
            # from the new center.
            nwafer = 19

            # Simulate the offset layout
            offset = qa.mult(
                qa.rotation(XAXIS, -np.radians(waferspace * platescale / 2)),
                qa.rotation(ZAXIS, -thirty),
            )
            wcenters = hex_layout(
                nwafer,
                4 * waferspace * platescale * u.degree,
                "",
                "",
                np.zeros(nwafer, dtype=np.float64) * u.degree,
                center=offset,
            )

            # Cut wafers outside the footprint
            centers = list()
            for p, q in wcenters.items():
                vec = qa.rotate(q["quat"], ZAXIS)
                dist = np.arccos(np.dot(vec, ZAXIS))
                if dist > np.radians(2.0 * waferspace * platescale):
                    continue
                centers.append(q["quat"])

            # Simulate each wafer, applying central rotation
            for windx, (wafer, center) in enumerate(zip(tubeprops["wafers"], centers)):
                wrot = qa.rotation(ZAXIS, wafer_ang_rad[windx])
                cent = qa.mult(center, wrot)
                dets = sim_wafer_detectors(
                    hw,
                    wafer,
                    platescale,
                    fwhm,
                    center=cent,
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
            np.zeros(91, dtype=np.float64) * u.degree,
            center=None,
        )

        for tindx, tube in enumerate(tubes):
            tubeprops = hw.data["tubes"][tube]
            waferspace = tubeprops["waferspace"]
            platescale = tubeprops["platescale"]
            location = f"{tubeprops['toast_hex_pos']:02}"
            wafer_ang_deg = tubeprops["wafer_angle"]
            wafer_ang_rad = np.radians(wafer_ang_deg)

            # NOTE:  If the CHLAT re-imaging optics are similar to SO LAT,
            # we may need a coordinate flip here for each tube.

            for windx, wafer in enumerate(tubeprops["wafers"]):
                wcenter = qa.mult(
                    qa.rotation(ZAXIS, wafer_ang_rad[windx]),
                    tcenters[location]["quat"],
                )
                dets = sim_wafer_detectors(
                    hw,
                    wafer,
                    platescale,
                    fwhm,
                    center=wcenter,
                )
                alldets.update(dets)

    if "detectors" in hw.data:
        hw.data["detectors"].update(alldets)
    else:
        hw.data["detectors"] = alldets
