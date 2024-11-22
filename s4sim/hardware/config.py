# Copyright (c) 2020-2020 CMB-S4 Collaboration.
# Full license can be found in the top level "LICENSE" file.
"""Hardware configuration utilities.
"""

import os
import re
import copy

from collections import OrderedDict

import gzip

import numpy as np

import toml


class Hardware(object):
    """Class representing a specific hardware configuration.

    The data is stored in a dictionary, and can be loaded / dumped to disk
    as well as trimmed to include only a subset of detectors.

    Args:
        path (str, optional): If specified, configuration is loaded from this
            file during construction.

    """

    def __init__(self, path=None):
        self.data = OrderedDict()
        if path is not None:
            self.load(path)

    def dump(self, path, overwrite=False, compress=False):
        """Write hardware config to a TOML file.

        Dump data to a TOML format file, optionally compressing the contents
        with gzip and optionally overwriting the file.

        Args:
            path (str): The file to write.
            overwrite (bool): If True, overwrite the file if it exists.
                If False, then existing files will cause an exception.
            compress (bool): If True, compress the data with gzip on write.

        Returns:
            None

        """
        if os.path.exists(path):
            if overwrite:
                os.remove(path)
            else:
                raise RuntimeError(
                    "Dump path {} already exists.  Use overwrite option".format(path)
                )
        if compress:
            with gzip.open(path, "wb") as f:
                dstr = toml.dumps(self.data)
                f.write(dstr.encode())
        else:
            with open(path, "w") as f:
                dstr = toml.dumps(self.data)
                f.write(dstr)
        return

    def load(self, path):
        """Read data from a TOML file.

        The file can either be regular text or a gzipped version of a TOML
        file.

        Args:
            path (str): The file to read.

        Returns:
            None

        """
        dstr = None
        try:
            with gzip.open(path, "rb") as f:
                dstr = f.read()
                self.data = toml.loads(dstr.decode())
        except OSError:
            with open(path, "r") as f:
                dstr = f.read()
                self.data = toml.loads(dstr)
        return

    def wafer_map(self):
        """Construct wafer mapping to other auxilliary data.

        Given the current data state, build dictionaries to go from wafers
        to all other non-detector info:  telescopes, tubes, cards, crates,
        and bands.  This is a convenient mapping when pruning the hardware
        information or doing other kinds of lookups.

        Returns:
            (dict): Nested dictionaries from wafers to other properties.

        """
        result = OrderedDict()

        tube_to_tele = dict()
        for tele, props in self.data["telescopes"].items():
            for tb in props["tubes"]:
                tube_to_tele[tb] = tele

        wafer_to_tube = dict()
        for tb, props in self.data["tubes"].items():
            for wf in props["wafers"]:
                wafer_to_tube[wf] = tb

        crate_to_card = dict()
        for crate, props in self.data["crates"].items():
            for card in props["cards"]:
                crate_to_card[card] = crate

        result["cards"] = {x: y["card"] for x, y in self.data["wafers"].items()}
        result["crates"] = {
            x: crate_to_card[y["card"]] for x, y in self.data["wafers"].items()
        }
        result["bands"] = {x: y["bands"] for x, y in self.data["wafers"].items()}
        result["tubes"] = wafer_to_tube
        result["telescopes"] = {
            x: tube_to_tele[wafer_to_tube[x]] for x in list(self.data["wafers"].keys())
        }
        return result

    def select(self, telescopes=None, tubes=None, match=dict()):
        """Select a subset of detectors.

        Select detectors whose properties match some criteria.  A new Hardware
        object is created and returned.  If a matching expression is not
        specified for a given property name, then this is equivalent to
        selecting all values of that property.

        Before selecting on detector properties, any telescope / tube filtering
        criteria are first applied.

        Each key of the "match" dictionary should be the name of a detector
        property to be considered for selection (e.g. band, wafer, pol, pixel).
        The value is a matching expression which can be:

            - A list of explicit values to match.
            - A string containing a regex expression to apply.

        Example:
            Imagine you wanted to select all 90GHz detectors on wafers 25 and
            26 which have "A" polarization and are located in pixels 20-29
            (recall the "." matches a single character)::

                new = hw.select(match={"wafer": ["25", "26"],
                                "band": "MF.1",
                                "pol": "A",
                                "pixel": "02."})

        Args:
            telescopes (str): A list of telescope names.
            tubes (str): A list of tube names.
            match (dict): The dictionary of property names and their matching
                expressions.

        Returns:
            (Hardware): A new Hardware instance with the selected detectors.

        """
        # First parse any telescope and tube options into a list of wafers
        wselect = None
        tbselect = None
        if telescopes is not None:
            tbselect = list()
            for tele in telescopes:
                tbselect.extend(self.data["telescopes"][tele]["tubes"])
        if tubes is not None:
            if tbselect is None:
                tbselect = list()
            tbselect.extend(tubes)
        if tbselect is not None:
            wselect = list()
            for tb in tbselect:
                wselect.extend(self.data["tubes"][tb]["wafers"])

        dets = self.data["detectors"]

        # Build regex matches for each property
        reg = dict()
        if "wafer" in match:
            # Handle wafer case separately, since we need to merge any
            # match with our telescope / tube selection of wafers above.
            k = "wafer"
            v = match[k]
            if wselect is None:
                # Just the regular behavior
                if isinstance(v, list):
                    reg[k] = re.compile(r"(^" + "$|^".join(v) + r"$)")
                else:
                    reg[k] = re.compile(v)
            else:
                # Merge our selection
                wall = list(wselect)
                if isinstance(v, list):
                    wall.extend(v)
                else:
                    wall.append(v)
                reg[k] = re.compile(r"(^" + "$|^".join(wall) + r"$)")
        elif wselect is not None:
            # No pattern in the match dictionary, just our list from the
            # telescope / tube selection.
            if wselect is None:
                reg["wafer"] = re.compile(r".*")
            else:
                reg["wafer"] = re.compile(r"(^" + "$|^".join(wselect) + r"$)")

        for k, v in match.items():
            if k == "wafer":
                # Already handled above
                continue
            else:
                if isinstance(v, list):
                    reg[k] = re.compile(r"(^" + "$|^".join(v) + r"$)")
                else:
                    reg[k] = re.compile(v)

        # Go through all detectors selecting things that match all fields
        newwafers = set()
        newdets = OrderedDict()
        for d, props in dets.items():
            keep = True
            for k, v in reg.items():
                if k in props:
                    test = v.match(props[k])
                    if test is None:
                        keep = False
                        break
            if keep:
                newwafers.add(props["wafer"])
                newdets[d] = copy.deepcopy(props)

        # Now compute the reduced set of auxilliary data needed for these
        # detectors.
        wafermap = self.wafer_map()

        # Copy this data
        hw = Hardware()
        hw.data = OrderedDict()
        for k, v in wafermap.items():
            hw.data[k] = OrderedDict()
            tocopy = set()
            for wf in newwafers:
                if isinstance(v[wf], list):
                    for iv in v[wf]:
                        tocopy.add(iv)
                else:
                    tocopy.add(v[wf])
            for elem in tocopy:
                hw.data[k][elem] = copy.deepcopy(self.data[k][elem])

        # Copy over the wafer data
        hw.data["wafers"] = OrderedDict()
        for wf in newwafers:
            hw.data["wafers"][wf] = copy.deepcopy(self.data["wafers"][wf])

        # And the detectors...
        hw.data["detectors"] = newdets

        return hw


def sim_nominal():
    """Return a simulated nominal hardware configuration.

    This returns a simulated Hardware object with the nominal instrument
    properties / metadata, but with an empty set of detector locations.
    This can then be passed to one of the detector simulation functions
    to build up the list of detectors.

    Returns:
        (Hardware): Hardware object with nominal metadata.
    """
    cnf = OrderedDict()

    bands = OrderedDict()

    bnd = OrderedDict()
    bnd["center"] = 20.0
    bnd["low"] = 17.5
    bnd["high"] = 22.5
    bnd["bandpass"] = ""
    bnd["NET"] = 406.0
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 1.0
    bnd["A"] = 0.07100
    bnd["C"] = 0.90757
    bnd["NET_corr"] = 1.28
    bnd["pwv_poly"] = 0.944529, 0.053562, 0.002317
    bands["CHLAT_f020"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 20.0
    bnd["low"] = 17.5
    bnd["high"] = 22.5
    bnd["bandpass"] = ""
    bnd["NET"] = 421.0
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 1.0
    bnd["A"] = 0.05443
    bnd["C"] = 0.92913
    bnd["NET_corr"] = 1.23
    bnd["pwv_poly"] = 0.989486, 0.032389, 0.001136
    bands["SPLAT_f020"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 24.75
    bnd["low"] = 21.5
    bnd["high"] = 28.0
    bnd["bandpass"] = ""
    bnd["NET"] = 384.0
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 1.0
    bnd["A"] = 0.06980
    bnd["C"] =  0.90909
    bnd["NET_corr"] = 1.30
    bnd["pwv_poly"] = 0.956524, 0.042312, 0.001480
    bands["CHLAT_f030"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 36.5
    bnd["low"] = 28.0
    bnd["high"] = 45.0
    bnd["bandpass"] = ""
    bnd["NET"] = 230.0
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 1.0
    bnd["A"] = 0.17386
    bnd["C"] = 0.77382
    bnd["NET_corr"] = 1.14
    bnd["pwv_poly"] = 0.987281, 0.012516, 0.000294
    bands["CHLAT_f040"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 24.75
    bnd["low"] = 21.5
    bnd["high"] = 28.0
    bnd["bandpass"] = ""
    bnd["NET"] = 407.0
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 1.0
    bnd["A"] = 0.06565
    bnd["C"] = 0.91454
    bnd["NET_corr"] = 1.27
    bnd["pwv_poly"] = 0.990859, 0.028229, 0.000774
    bands["SPLAT_f030"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 36.5
    bnd["low"] = 28.0
    bnd["high"] = 45.0
    bnd["bandpass"] = ""
    bnd["NET"] = 260.0
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 1.0
    bnd["A"] = 0.26014
    bnd["C"] = 0.66203
    bnd["NET_corr"] = 1.17
    bnd["pwv_poly"] = 0.995952, 0.012515, 0.000299
    bands["SPLAT_f040"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 24.75
    bnd["low"] = 21.5
    bnd["high"] = 28.0
    bnd["bandpass"] = ""
    bnd["NET"] = 214.0
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 1.0
    bnd["A"] = 0.13706
    bnd["C"] = 0.82200
    bnd["NET_corr"] = 1.06
    bnd["pwv_poly"] = 0.920270, 0.076947, 0.003369
    bands["SAT_f030"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 36.5
    bnd["low"] = 28.0
    bnd["high"] = 45.0
    bnd["bandpass"] = ""
    bnd["NET"] = 148.0
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 1.0
    bnd["A"] = 0.23230
    bnd["C"] = 0.69828
    bnd["NET_corr"] = 1.04
    bnd["pwv_poly"] = 0.983677, 0.016068, 0.000373
    bands["SAT_f040"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 91.5
    bnd["low"] = 77.0
    bnd["high"] = 106.0
    bnd["bandpass"] = ""
    bnd["NET"] = 288.0
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 1.0
    bnd["A"] = 0.17224
    bnd["C"] = 0.77562
    bnd["NET_corr"] = 1.18
    bnd["pwv_poly"] = 0.950153, 0.048503, 0.001707
    bands["CHLAT_f090"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 148.5
    bnd["low"] = 128.0
    bnd["high"] = 169.0
    bnd["bandpass"] = ""
    bnd["NET"] = 313.0
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 1.0
    bnd["A"] = 0.21507
    bnd["C"] = 0.71969
    bnd["NET_corr"] = 1.01
    bnd["pwv_poly"] = 0.804635, 0.188923, 0.007874
    bands["CHLAT_f150"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 91.5
    bnd["low"] = 77.0
    bnd["high"] = 106.0
    bnd["bandpass"] = ""
    bnd["NET"] = 324.0
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 1.0
    bnd["A"] = 0.23242
    bnd["C"] = 0.69750
    bnd["NET_corr"] = 1.18
    bnd["pwv_poly"] = 0.983165, 0.051852, 0.001852
    bands["SPLAT_f090"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 148.5
    bnd["low"] = 128.0
    bnd["high"] = 169.0
    bnd["bandpass"] = ""
    bnd["NET"] = 297.0
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 1.0
    bnd["A"] = 0.19995
    bnd["C"] = 0.73946
    bnd["NET_corr"] = 1.01
    bnd["pwv_poly"] = 0.921211, 0.241889, 0.011086
    bands["SPLAT_f150"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 85.0
    bnd["low"] = 74.8
    bnd["high"] = 95.2
    bnd["bandpass"] = ""
    bnd["NET"] = 245.0
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 1.0
    bnd["A"] = 0.27782
    bnd["C"] = 0.63888
    bnd["NET_corr"] = 1.04
    bnd["pwv_poly"] = 0.953624, 0.045185, 0.001529
    bands["SAT_f085"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 145.0
    bnd["low"] = 129.1
    bnd["high"] = 161.0
    bnd["bandpass"] = ""
    bnd["NET"] = 302.0
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 1.0
    bnd["A"] = 0.21203
    bnd["C"] = 0.72385
    bnd["NET_corr"] = 1.01
    bnd["pwv_poly"] = 0.845158, 0.150156, 0.005818
    bands["SAT_f145"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 95.0
    bnd["low"] = 83.6
    bnd["high"] = 106.4
    bnd["bandpass"] = ""
    bnd["NET"] = 228.0
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 1.0
    bnd["A"] = 0.22131
    bnd["C"] = 0.71200
    bnd["NET_corr"] = 1.03
    bnd["pwv_poly"] = 0.938631, 0.059726, 0.002090
    bands["SAT_f095"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 155.0
    bnd["low"] = 138.0
    bnd["high"] = 172.1
    bnd["bandpass"] = ""
    bnd["NET"] = 345.0
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 1.0
    bnd["A"] = 0.23468
    bnd["C"] = 0.69439
    bnd["NET_corr"] = 1.00
    bnd["pwv_poly"] = 0.792946, 0.200998, 0.007569
    bands["SAT_f155"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 91.5
    bnd["low"] = 77.0
    bnd["high"] = 106.0
    bnd["bandpass"] = ""
    bnd["NET"] = 221.0
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 1.0
    bnd["A"] = 0.22609
    bnd["C"] = 0.70593
    bnd["NET_corr"] = 1.03
    bnd["pwv_poly"] = 0.947923, 0.050705, 0.001751
    bands["SAT_f090"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 148.5
    bnd["low"] = 128.0
    bnd["high"] = 169.0
    bnd["bandpass"] = ""
    bnd["NET"] = 287.0
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 1.0
    bnd["A"] = 0.22482
    bnd["C"] = 0.70722
    bnd["NET_corr"] = 1.03
    bnd["pwv_poly"] = 0.826483, 0.168378, 0.006407
    bands["SAT_f150"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 227.0
    bnd["low"] = 198.0
    bnd["high"] = 256.0
    bnd["bandpass"] = ""
    bnd["NET"] = 670.0
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 1.0
    bnd["A"] = 0.32681
    bnd["C"] = 0.57429
    bnd["NET_corr"] = 1.02
    bnd["pwv_poly"] = 0.675399, 0.307360, 0.019666
    bands["CHLAT_f220"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 285.5
    bnd["low"] = 256.0
    bnd["high"] = 315.0
    bnd["bandpass"] = ""
    bnd["NET"] = 1615.0
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 1.0
    bnd["A"] = 0.40724
    bnd["C"] =  0.46986
    bnd["NET_corr"] = 1.01
    bnd["pwv_poly"] = 0.582919, 0.387797, 0.032451
    bands["CHLAT_f280"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 227.0
    bnd["low"] = 198.0
    bnd["high"] = 256.0
    bnd["bandpass"] = ""
    bnd["NET"] = 580.0
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 1.0
    bnd["A"] = 0.27310
    bnd["C"] = 0.64416
    bnd["NET_corr"] = 1.02
    bnd["pwv_poly"] = 0.854996, 0.441083, 0.033158
    bands["SPLAT_f220"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 285.5
    bnd["low"] = 256.0
    bnd["high"] = 315.0
    bnd["bandpass"] = ""
    bnd["NET"] = 1327.0
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 1.0
    bnd["A"] = 0.33182
    bnd["C"] = 0.56786
    bnd["NET_corr"] = 1.00
    bnd["pwv_poly"] = 0.802794, 0.594632, 0.061426
    bands["SPLAT_f280"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 227.0
    bnd["low"] = 198.0
    bnd["high"] = 256.0
    bnd["bandpass"] = ""
    bnd["NET"] = 720.0
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 1.0
    bnd["A"] = 0.30314
    bnd["C"] = 0.60552
    bnd["NET_corr"] = 1.01
    bnd["pwv_poly"] = 0.728284, 0.258698, 0.015039
    bands["SAT_f220"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 285.5
    bnd["low"] = 256.0
    bnd["high"] = 315.0
    bnd["bandpass"] = ""
    bnd["NET"] = 1817.0
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 1.0
    bnd["A"] = 0.36697
    bnd["C"] = 0.52293
    bnd["NET_corr"] = 1.00
    bnd["pwv_poly"] = 0.658761, 0.319594, 0.024221
    bands["SAT_f280"] = bnd

    cnf["bands"] = bands

    wafers = OrderedDict()

    wtypes = ["CHLAT_ULF","SPLAT_ULF", "CHLAT_LF", "SPLAT_LF", "SAT_LF", "CHLAT_MF", "SPLAT_MF", "SAT_MFL", "SAT_MFH", "SAT_MF", "CHLAT_HF", "SPLAT_HF", "SAT_HF"]
    wcnt = {
        "CHLAT_ULF": 4*2,
        "SPLAT_ULF": 4,
        "CHLAT_LF": 8*2,
        "SPLAT_LF": 9,
        "SAT_LF": (12) * 1,
        "CHLAT_MF": 54*2,
        "SPLAT_MF": 54,
        "SAT_MFL": 0,
        "SAT_MFH": 0,
        "SAT_MF": (12) * 6,
        "CHLAT_HF": 19*2,
        "SPLAT_HF": 18,
        "SAT_HF": (12) * 2,
    }
    wpac = {
    	"CHLAT_ULF": "RP",
        "SPLAT_ULF": "RP",
        "CHLAT_LF": "RP",
        "SPLAT_LF": "RP",
        "SAT_LF": "HP",
        "CHLAT_MF": "RP",
        "SPLAT_MF": "RP",
        "SAT_MFL": "HP",
        "SAT_MFH": "HP",
        "SAT_MF": "RP",
        "CHLAT_HF": "HP",
        "SPLAT_HF": "HP",
        "SAT_HF": "HP",
    }
    wnp = {
        "CHLAT_ULF": 27,
        "SPLAT_ULF": 27,
        "CHLAT_LF": 48,
        "SPLAT_LF": 48,
        "SAT_LF": 37,
        "CHLAT_MF": 432,
        "SPLAT_MF": 432,
        "SAT_MFL": 217,
        "SAT_MFH": 271,
        "SAT_MF": 432,
        "CHLAT_HF": 469,
        "SPLAT_HF": 469,
        "SAT_HF": 469,
    }
    wpixmm = {
        "CHLAT_ULF": 21.1,
        "SPLAT_ULF": 21.1,
        "CHLAT_LF": 16.1,
        "SPLAT_LF": 16.1,
        "SAT_LF": 19.1,
        "CHLAT_MF": 5.3,
        "SPLAT_MF": 5.3,
        "SAT_MFL": 7.65,
        "SAT_MFH": 6.85,
        "SAT_MF": 5.3,
        "CHLAT_HF": 5.2,
        "SPLAT_HF": 5.2,
        "SAT_HF": 5.2,
    }
    wrhombgap = {
        "CHLAT_ULF": 2.827,
        "SPLAT_ULF": 2.827,
        "CHLAT_LF": 2.157,
        "SPLAT_LF": 2.157,
        "SAT_LF": 4.167,
        "CHLAT_MF": 0.71,
        "SPLAT_MF": 0.71,
        "SAT_MFL": 1.273,
        "SAT_MFH": 0.71,
        "SAT_MF": 0.71,
        "CHLAT_HF": 0.71,
        "SPLAT_HF": 0.71,
        "SAT_HF": 0.71,
    }
    wbd = {
        "CHLAT_ULF": ["CHLAT_f020"],
        "SPLAT_ULF": ["SPLAT_f020"],
        "CHLAT_LF": ["CHLAT_f030", "CHLAT_f040"],
        "SPLAT_LF": ["SPLAT_f030", "SPLAT_f040"],
        "SAT_LF": ["SAT_f030", "SAT_f040"],
        "CHLAT_MF": ["CHLAT_f090", "CHLAT_f150"],
        "SPLAT_MF": ["SPLAT_f090", "SPLAT_f150"],
        "SAT_MFL": ["SAT_f085", "SAT_f145"],
        "SAT_MFH": ["SAT_f095", "SAT_f155"],
        "SAT_MF": ["SAT_f090", "SAT_f150"],
        "CHLAT_HF": ["CHLAT_f220", "CHLAT_f280"],
        "SPLAT_HF": ["SPLAT_f220", "SPLAT_f280"],
        "SAT_HF": ["SAT_f220", "SAT_f280"],
    }
    #location positions of mechanical pin holes that kill pixels on the wafers
    pins = {
        "CHLAT_ULF": [],
        "SPLAT_ULF": [],
        "CHLAT_LF": [],
        "SPLAT_LF": [],
        "SAT_LF": [],
        "CHLAT_MF": [210,220],
        "SPLAT_MF": [210,220],
        "SAT_MFL": [0,127],
        "SAT_MFH": [0,169],
        "SAT_MF": [210,220],
        "CHLAT_HF": [0,331],
        "SPLAT_HF": [0,331],
        "SAT_HF": [0,331],
    }
    windx = 0
    cardindx = 0
    for wt in wtypes:
        for ct in range(wcnt[wt]):
            wn = "{:03d}".format(windx)
            wf = OrderedDict()
            wf["type"] = wt
            wf["packing"] = wpac[wt]
            wf["rhombusgap"] = wrhombgap[wt]
            wf["npixel"] = wnp[wt]
            wf["pixsize"] = wpixmm[wt]
            wf["bands"] = wbd[wt]
            wf["card"] = "{:02d}".format(cardindx)
            wf["pins"] = pins[wt]
            cardindx += 1
            wafers[wn] = wf
            windx += 1

    cnf["wafers"] = wafers

    tubes = OrderedDict()

    woff = {
        "CHLAT_ULF": 0,
        "SPLAT_ULF": 0,
        "CHLAT_LF": 0,
        "SPLAT_LF": 0,
        "SAT_LF": 0,
        "CHLAT_MF": 0,
        "SPLAT_MF": 0,
        "SAT_MFL": 0,
        "SAT_MFH": 0,
        "SAT_MF": 0,
        "CHLAT_HF": 0,
        "SPLAT_HF": 0,
        "SAT_HF": 0,
    }

    # added in tube platescale because SAT HF has a different platescale
    # Wafers are arranges in the tube with the tube platescale
    # Telescope platescale is used for spacing tubes
    ltubes = [
        "CHLAT_HF",
        "CHLAT_HF",
        "CHLAT_HF",
        "CHLAT_HF",
        "CHLAT_HF",
        "CHLAT_HF",
        "CHLAT_HF",
        "CHLAT_HF",
        "CHLAT_HF",
        "CHLAT_HF",
        "CHLAT_HF",
        "CHLAT_HF",
        "CHLAT_HF",
        "CHLAT_HF",
        "CHLAT_HF",
        "CHLAT_HF",
        "CHLAT_HF",
        "CHLAT_HF",
        "CHLAT_HF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_LF",
        "CHLAT_MF",
        "CHLAT_LF",
        "CHLAT_ULF",
        "CHLAT_LF",
        "CHLAT_ULF",
        "CHLAT_LF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_LF",
        "CHLAT_MF",
        "CHLAT_LF",
        "CHLAT_ULF",
        "CHLAT_LF",
        "CHLAT_ULF",
        "CHLAT_LF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_HF",
        "CHLAT_HF",
        "CHLAT_HF",
        "CHLAT_HF",
        "CHLAT_HF",
        "CHLAT_HF",
        "CHLAT_HF",
        "CHLAT_HF",
        "CHLAT_HF",
        "CHLAT_HF",
        "CHLAT_HF",
        "CHLAT_HF",
        "CHLAT_HF",
        "CHLAT_HF",
        "CHLAT_HF",
        "CHLAT_HF",
        "CHLAT_HF",
        "CHLAT_HF",
        "CHLAT_HF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_LF",
        "CHLAT_MF",
        "CHLAT_LF",
        "CHLAT_ULF",
        "CHLAT_LF",
        "CHLAT_ULF",
        "CHLAT_LF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_LF",
        "CHLAT_MF",
        "CHLAT_LF",
        "CHLAT_ULF",
        "CHLAT_LF",
        "CHLAT_ULF",
        "CHLAT_LF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "SPLAT_MF",
        "SPLAT_MF",
        "SPLAT_MF",
        "SPLAT_MF",
        "SPLAT_MF",
        "SPLAT_MF",
        "SPLAT_MF",
        "SPLAT_HF",
        "SPLAT_ULF",
        "SPLAT_HF",
        "SPLAT_MF",
        "SPLAT_HF",
        "SPLAT_ULF",
        "SPLAT_HF",
        "SPLAT_ULF",
        "SPLAT_HF",
        "SPLAT_MF",
        "SPLAT_HF",
        "SPLAT_ULF",
        "SPLAT_MF",
        "SPLAT_MF",
        "SPLAT_MF",
        "SPLAT_MF",
        "SPLAT_MF",
        "SPLAT_MF",
        "SPLAT_MF",
        "SPLAT_MF",
        "SPLAT_MF",
        "SPLAT_MF",
        "SPLAT_MF",
        "SPLAT_MF",
        "SPLAT_MF",
        "SPLAT_MF",
        "SPLAT_MF",
        "SPLAT_MF",
        "SPLAT_MF",
        "SPLAT_MF",
        "SPLAT_HF",
        "SPLAT_MF",
        "SPLAT_HF",
        "SPLAT_MF",
        "SPLAT_HF",
        "SPLAT_MF",
        "SPLAT_HF",
        "SPLAT_MF",
        "SPLAT_HF",
        "SPLAT_MF",
        "SPLAT_HF",
        "SPLAT_MF",
        "SPLAT_HF",
        "SPLAT_MF",
        "SPLAT_HF",
        "SPLAT_MF",
        "SPLAT_HF",
        "SPLAT_MF",
        "SPLAT_HF",
        "SPLAT_MF",
        "SPLAT_HF",
        "SPLAT_MF",
        "SPLAT_HF",
        "SPLAT_MF",
        "SPLAT_MF",
        "SPLAT_LF",
        "SPLAT_MF",
        "SPLAT_LF",
        "SPLAT_MF",
        "SPLAT_LF",
        "SPLAT_MF",
        "SPLAT_MF",
        "SPLAT_LF",
        "SPLAT_MF",
        "SPLAT_MF",
        "SPLAT_LF",
        "SPLAT_MF",
        "SPLAT_MF",
        "SPLAT_LF",
        "SPLAT_MF",
        "SPLAT_MF",
        "SPLAT_LF",
        "SPLAT_MF",
        "SPLAT_MF",
        "SPLAT_LF",
        "SPLAT_MF",
        "SPLAT_MF",
        "SPLAT_LF",
    ]
    # TOAST hexagon layout positions in Xi/Eta coordinates
    ltube_toasthex_pos = [
        0,
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        21,
        22,
        23,
        24,
        25,
        26,
        27,
        28,
        29,
        30,
        31,
        32,
        33,
        34,
        35,
        36,
        37,
        38,
        39,
        40,
        41,
        42,
        43,
        44,
        45,
        46,
        47,
        48,
        49,
        50,
        51,
        52,
        53,
        54,
        55,
        56,
        57,
        58,
        59,
        60,
        62,
        63,
        64,
        65,
        67,
        68,
        69,
        70,
        72,
        73,
        74,
        75,
        77,
        78,
        79,
        80,
        82,
        83,
        84,
        85,
        87,
        88,
        89,
        90,
        0,
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        21,
        22,
        23,
        24,
        25,
        26,
        27,
        28,
        29,
        30,
        31,
        32,
        33,
        34,
        35,
        36,
        37,
        38,
        39,
        40,
        41,
        42,
        43,
        44,
        45,
        46,
        47,
        48,
        49,
        50,
        51,
        52,
        53,
        54,
        55,
        56,
        57,
        58,
        59,
        60,
        62,
        63,
        64,
        65,
        67,
        68,
        69,
        70,
        72,
        73,
        74,
        75,
        77,
        78,
        79,
        80,
        82,
        83,
        84,
        85,
        87,
        88,
        89,
        90,
        0,
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        21,
        22,
        23,
        24,
        25,
        26,
        27,
        28,
        29,
        30,
        31,
        32,
        33,
        34,
        35,
        36,
        37,
        38,
        39,
        40,
        41,
        42,
        43,
        44,
        45,
        46,
        47,
        48,
        49,
        50,
        51,
        52,
        53,
        54,
        55,
        56,
        57,
        58,
        59,
        60,
        62,
        63,
        64,
        65,
        67,
        68,
        69,
        70,
        72,
        73,
        74,
        75,
        77,
        78,
        79,
        80,
        82,
        83,
        84,
        85,
        87,
        88,
        89,
        90
    ]
    # tindx is the tube number we have 85*3=255
    for tindx in range(255):
        nm = "LT{:d}".format(tindx)
        ttyp = ltubes[tindx]
        tb = OrderedDict()
        tb["type"] = ttyp
        tb["waferspace"] = 124.
        tb["wafers"] = list()
        # There is currently one wafer per tube.
        for tw in range(1):
            off = 0
            for w, props in cnf["wafers"].items():
                if props["type"] == ttyp:
                    if off == woff[ttyp]:
                        tb["wafers"].append(w)
                        woff[ttyp] += 1
                        break
                    off += 1
        if ttyp == "CHLAT_HF" or ttyp == "SPLAT_HF":
            # These are hex wafers
            tb["wafer_angle"] = [-90.0 for tw in range(1)] # Degrees
        else:
            # These are rhombi-hex wafers
            tb["wafer_angle"] = [0.0 for tw in range(1)] # Degrees
        if tindx < 170:
            # CHLAT platescale
            tb["platescale"] = 0.003964
        else:
            # SPLAT platescale
            tb["platescale"] = 0.00429
        tb["toast_hex_pos"] = ltube_toasthex_pos[tindx]
        tubes[nm] = tb

    stubes = [
        "SAT_MF",
        "SAT_MF",
        "SAT_HF",
        "SAT_MF",
        "SAT_MF",
        "SAT_HF",
        "SAT_MF",
        "SAT_MF",
        "SAT_LF",
    ]
    stube_toasthex_pos = [0, 0, 0, 0, 0, 0, 0, 0, 0]
    for tindx in range(9):
        nm = "ST{:d}".format(tindx)
        ttyp = stubes[tindx]
        tb = OrderedDict()
        tb["type"] = ttyp
        tb["waferspace"] = 122.0
        tb["wafers"] = list()
        if ttyp == "SAT_HF":
            for tw in range(12):
                off = 0
                for w, props in cnf["wafers"].items():
                    if props["type"] == ttyp:
                        if off == woff[ttyp]:
                            tb["wafers"].append(w)
                            woff[ttyp] += 1
                            break
                        off += 1
            tb["wafer_angle"] = [ # Degrees
                30.0,
                30.0,
                30.0,
                -150.0,
                -150.0,
                30.0,
                30.0,
                30.0,
                -30.0,
                90.0,
                30.0,
                30.0,
            ]
            # 30 deg, 10008 detectors/band
            tb["platescale"] = 0.070093/0.9909
            tb["FOV_cut"] = 34.9
        elif ttyp == "SAT_MFL":
            for tw in range(12):
                off = 0
                for w, props in cnf["wafers"].items():
                    if props["type"] == ttyp:
                        if off == woff[ttyp]:
                            tb["wafers"].append(w)
                            woff[ttyp] += 1
                            break
                        off += 1
            tb["wafer_angle"] = [ # Degrees
                30.0,
                30.0,
                30.0,
                -150.0,
                -150.0,
                30.0,
                30.0,
                30.0,
                -30.0,
                90.0,
                30.0,
                30.0,
            ]
            # 30 deg, 3048 detectors/band
            tb["platescale"] = 0.070093/0.9905
            tb["FOV_cut"] = 30.0
        elif ttyp == "SAT_MFH":
            for tw in range(12):
                off = 0
                for w, props in cnf["wafers"].items():
                    if props["type"] == ttyp:
                        if off == woff[ttyp]:
                            tb["wafers"].append(w)
                            woff[ttyp] += 1
                            break
                        off += 1
            tb["wafer_angle"] = [ # Degrees
                30.0,
                30.0,
                30.0,
                -150.0,
                -150.0,
                30.0,
                30.0,
                30.0,
                -30.0,
                90.0,
                30.0,
                30.0,
            ]
            # 30 deg, 3552 det/band
            tb["platescale"] = 0.070093/0.9931
            tb["FOV_cut"] = 30.0
        elif ttyp == "SAT_MF":
            for tw in range(12):
                off = 0
                for w, props in cnf["wafers"].items():
                    if props["type"] == ttyp:
                        if off == woff[ttyp]:
                            tb["wafers"].append(w)
                            woff[ttyp] += 1
                            break
                        off += 1
            tb["wafer_angle"] = [ # Degrees
                0.0,
                0.0,
                0.0,
                -60.0,
                -60.0,
                0.0,
                0.0,
                60.0,
                60.0,
                180.0,
                180.0,
                0.0,
            ]
            # 30 deg, 3552 det/band
            tb["platescale"] = 0.070093/0.9931
            tb["FOV_cut"] = 34.9
        else:
            for tw in range(12):
                off = 0
                for w, props in cnf["wafers"].items():
                    if props["type"] == ttyp:
                        if off == woff[ttyp]:
                            tb["wafers"].append(w)
                            woff[ttyp] += 1
                            break
                        off += 1
            tb["wafer_angle"] = [ # Degrees
                30.0,
                30.0,
                30.0,
                -150.0,
                -150.0,
                30.0,
                30.0,
                30.0,
                -30.0,
                90.0,
                30.0,
                30.0,
            ]
            # 30 deg, 252 det/band
            tb["platescale"] = 0.070093/0.9855
            tb["FOV_cut"] = 34.9
        tb["toast_hex_pos"] = stube_toasthex_pos[tindx]
        tubes[nm] = tb

    cnf["tubes"] = tubes

    telescopes = OrderedDict()

    lfwhm = OrderedDict()
    lfwhm["CHLAT_f020"] = 9.6
    lfwhm["SPLAT_f020"] = 11.2
    lfwhm["CHLAT_f030"] = 7.8
    lfwhm["SPLAT_f030"] = 9.1
    lfwhm["CHLAT_f040"] = 5.3
    lfwhm["SPLAT_f040"] = 6.2
    lfwhm["CHLAT_f090"] = 2.1
    lfwhm["SPLAT_f090"] = 2.5
    lfwhm["CHLAT_f150"] = 1.3
    lfwhm["SPLAT_f150"] = 1.5
    lfwhm["CHLAT_f220"] = 0.95
    lfwhm["SPLAT_f220"] = 1.1
    lfwhm["CHLAT_f280"] = 0.83
    lfwhm["SPLAT_f280"] = 1.0

    tele = OrderedDict()
    tele["tubes"] = [
        "LT0",
        "LT1",
        "LT2",
        "LT3",
        "LT4",
        "LT5",
        "LT6",
        "LT7",
        "LT8",
        "LT9",
        "LT10",
        "LT11",
        "LT12",
        "LT13",
        "LT14",
        "LT15",
        "LT16",
        "LT17",
        "LT18",
        "LT19",
        "LT20",
        "LT21",
        "LT22",
        "LT23",
        "LT24",
        "LT25",
        "LT26",
        "LT27",
        "LT28",
        "LT29",
        "LT30",
        "LT31",
        "LT32",
        "LT33",
        "LT34",
        "LT35",
        "LT36",
        "LT37",
        "LT38",
        "LT39",
        "LT40",
        "LT41",
        "LT42",
        "LT43",
        "LT44",
        "LT45",
        "LT46",
        "LT47",
        "LT48",
        "LT49",
        "LT50",
        "LT51",
        "LT52",
        "LT53",
        "LT54",
        "LT55",
        "LT56",
        "LT57",
        "LT58",
        "LT59",
        "LT60",
        "LT61",
        "LT62",
        "LT63",
        "LT64",
        "LT65",
        "LT66",
        "LT67",
        "LT68",
        "LT69",
        "LT70",
        "LT71",
        "LT72",
        "LT73",
        "LT74",
        "LT75",
        "LT76",
        "LT77",
        "LT78",
        "LT79",
        "LT80",
        "LT81",
        "LT82",
        "LT83",
        "LT84"
    ]
    tele["platescale"] = 0.003964
    # This tube spacing in mm corresponds to 0.87 degrees projected on
    # the sky at a plate scale of 0.003964 deg/mm
    # The physical tube spacing may be different.
    tele["tubespace"] =  0.87 / 0.003964
    tele["fwhm"] = lfwhm
    telescopes["LAT0"] = tele

    tele = OrderedDict()
    tele["tubes"] = [
        "LT85",
        "LT86",
        "LT87",
        "LT88",
        "LT89",
        "LT90",
        "LT91",
        "LT92",
        "LT93",
        "LT94",
        "LT95",
        "LT96",
        "LT97",
        "LT98",
        "LT99",
        "LT100",
        "LT101",
        "LT102",
        "LT103",
        "LT104",
        "LT105",
        "LT106",
        "LT107",
        "LT108",
        "LT109",
        "LT110",
        "LT111",
        "LT112",
        "LT113",
        "LT114",
        "LT115",
        "LT116",
        "LT117",
        "LT118",
        "LT119",
        "LT120",
        "LT121",
        "LT122",
        "LT123",
        "LT124",
        "LT125",
        "LT126",
        "LT127",
        "LT128",
        "LT129",
        "LT130",
        "LT131",
        "LT132",
        "LT133",
        "LT134",
        "LT135",
        "LT136",
        "LT137",
        "LT138",
        "LT139",
        "LT140",
        "LT141",
        "LT142",
        "LT143",
        "LT144",
        "LT145",
        "LT146",
        "LT147",
        "LT148",
        "LT149",
        "LT150",
        "LT151",
        "LT152",
        "LT153",
        "LT154",
        "LT155",
        "LT156",
        "LT157",
        "LT158",
        "LT159",
        "LT160",
        "LT161",
        "LT162",
        "LT163",
        "LT164",
        "LT165",
        "LT166",
        "LT167",
        "LT168",
        "LT169"
    ]
    tele["platescale"] =  0.003964
    # This tube spacing in mm corresponds to 0.87 degrees projected on
    # the sky at a plate scale of 0.003964 deg/mm
    # The physical tube spacing may be different.
    tele["tubespace"] =  0.87 / 0.003964
    tele["fwhm"] = lfwhm
    telescopes["LAT1"] = tele

    tele = OrderedDict()
    tele["tubes"] = [
        "LT170",
        "LT171",
        "LT172",
        "LT173",
        "LT174",
        "LT175",
        "LT176",
        "LT177",
        "LT178",
        "LT179",
        "LT180",
        "LT181",
        "LT182",
        "LT183",
        "LT184",
        "LT185",
        "LT186",
        "LT187",
        "LT188",
        "LT189",
        "LT190",
        "LT191",
        "LT192",
        "LT193",
        "LT194",
        "LT195",
        "LT196",
        "LT197",
        "LT198",
        "LT199",
        "LT200",
        "LT201",
        "LT202",
        "LT203",
        "LT204",
        "LT205",
        "LT206",
        "LT207",
        "LT208",
        "LT209",
        "LT210",
        "LT211",
        "LT212",
        "LT213",
        "LT214",
        "LT215",
        "LT216",
        "LT217",
        "LT218",
        "LT219",
        "LT220",
        "LT221",
        "LT222",
        "LT223",
        "LT224",
        "LT225",
        "LT226",
        "LT227",
        "LT228",
        "LT229",
        "LT230",
        "LT231",
        "LT232",
        "LT233",
        "LT234",
        "LT235",
        "LT236",
        "LT237",
        "LT238",
        "LT239",
        "LT240",
        "LT241",
        "LT242",
        "LT243",
        "LT244",
        "LT245",
        "LT246",
        "LT247",
        "LT248",
        "LT249",
        "LT250",
        "LT251",
        "LT252",
        "LT253",
        "LT254"
    ]
    tele["platescale"] = 0.00429
    # This tube spacing in mm corresponds to 0.94 degrees projected on
    # the sky at a plate scale of 0.00429 deg/mm
    # The physical tube spacing may be different.
    tele["tubespace"] =  0.94 / 0.00429
    tele["fwhm"] = lfwhm
    telescopes["LAT2"] = tele

    #SAT beams
    sfwhm = OrderedDict()
    sfwhm["SAT_f030"] = 79.2
    sfwhm["SAT_f040"] = 59.4
    sfwhm["SAT_f085"] = 23.6
    sfwhm["SAT_f145"] = 15.0
    sfwhm["SAT_f095"] = 21.2
    sfwhm["SAT_f155"] = 13.9
    sfwhm["SAT_f090"] = 21.4
    sfwhm["SAT_f150"] = 14.0
    sfwhm["SAT_f220"] = 9.4
    sfwhm["SAT_f280"] = 8.3

    tele = OrderedDict()
    tele["tubes"] = ["ST0", "ST1", "ST2"]
    tele["platescale"] = 0.056689
    tele["tubespace"] = 700.0
    tele["fwhm"] = sfwhm
    telescopes["SAT1"] = tele

    tele = OrderedDict()
    tele["tubes"] = ["ST3", "ST4", "ST5"]
    tele["platescale"] = 0.056689
    tele["tubespace"] = 700.0
    tele["fwhm"] = sfwhm
    telescopes["SAT2"] = tele

    tele = OrderedDict()
    tele["tubes"] = ["ST6", "ST7", "ST8"]
    tele["platescale"] = 0.056689
    tele["tubespace"] = 700.0
    tele["fwhm"] = sfwhm
    telescopes["SAT3"] = tele

    cnf["telescopes"] = telescopes

    cards = OrderedDict()
    crates = OrderedDict()

    crt_indx = 0

    for tel in cnf["telescopes"]:
        crn = "{:d}".format(crt_indx)
        crt = OrderedDict()
        crt["cards"] = list()
        crt["telescope"] = tel

        ## get all the wafer card numbers for a telescope
        tb_wfrs = [cnf["tubes"][t]["wafers"] for t in cnf["telescopes"][tel]["tubes"]]
        tl_wfrs = [i for sl in tb_wfrs for i in sl]
        wafer_cards = [cnf["wafers"][w]["card"] for w in tl_wfrs]

        # add all cards to the card table and assign to crates
        for crd in wafer_cards:
            cdprops = OrderedDict()
            cdprops["nbias"] = 12
            cdprops["ncoax"] = 2
            cdprops["nchannel"] = 2000
            cards[crd] = cdprops

            crt["cards"].append(crd)

            # name new crates when current one is full
            # 6 cards/crate, in future for partial wafers change card number accordingly
            if len(crt["cards"]) >= 6:
                crates[crn] = crt
                crt_indx += 1
                crn = "{:d}".format(crt_indx)
                crt = OrderedDict()
                crt["cards"] = list()
                crt["telescope"] = tel

        # each telescope starts with a new crate
        crates[crn] = crt
        crt_indx += 1

    cnf["cards"] = cards
    cnf["crates"] = crates

    # Add an empty set of detectors here, in case the user just wants access to
    # the hardware metadata.
    cnf["detectors"] = OrderedDict()

    hw = Hardware()
    hw.data = cnf

    return hw
