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
            telescopes (str): A regex string to apply to telescope names or a
                list of explicit names.
            tubes (str): A regex string to apply to tube names or a list of
                explicit names.
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


def get_example():
    """Return an example Hardware config with the required sections.

    The returned Hardware object has 4 fake detectors as an example.  These
    detectors can be replaced by the results of other simulation functions.

    Returns:
        (Hardware): Hardware object with example parameters.

    """
    cnf = OrderedDict()

    bands = OrderedDict()

    bnd = OrderedDict()
    bnd["center"] = 20.0
    bnd["low"] = 17.5
    bnd["high"] = 22.5
    bnd["bandpass"] = ""
    bnd["NET"] = 332.7
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 3.5
    bnd["A"] = 0.08386
    bnd["C"] = 0.89063
    bnd["NET_corr"] = 1.22
    bnd["pwv_poly"] = 0.971399, 0.088721, 0.001181
    bands["SPLAT_f020"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 25.75
    bnd["low"] = 21.5
    bnd["high"] = 30.0
    bnd["bandpass"] = ""
    bnd["NET"] = 307.7
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 3.5
    bnd["A"] = 0.08426
    bnd["C"] = 0.89010
    bnd["NET_corr"] = 1.28
    bnd["pwv_poly"] = 0.934479, 0.064763, 0.001228
    bands["CHLAT_f030"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 38.75
    bnd["low"] = 30.0
    bnd["high"] = 47.5
    bnd["bandpass"] = ""
    bnd["NET"] = 240.7
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 3.5
    bnd["A"] = 0.24749
    bnd["C"] = 0.67744
    bnd["NET_corr"] = 1.04
    bnd["pwv_poly"] = 0.982523, 0.017230, 0.000373
    bands["CHLAT_f040"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 25.75
    bnd["low"] = 21.5
    bnd["high"] = 30.0
    bnd["bandpass"] = ""
    bnd["NET"] = 286.7
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 3.5
    bnd["A"] = 0.10372
    bnd["C"] = 0.86476
    bnd["NET_corr"] = 1.24
    bnd["pwv_poly"] = 0.981351, 0.057846, 0.000779
    bands["SPLAT_f030"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 38.75
    bnd["low"] = 30.0
    bnd["high"] = 47.5
    bnd["bandpass"] = ""
    bnd["NET"] = 269.0
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 3.5
    bnd["A"] = 0.38043
    bnd["C"] = 0.50426
    bnd["NET_corr"] = 1.04
    bnd["pwv_poly"] = 0.994233, 0.017840, 0.000388
    bands["SPLAT_f040"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 25.75
    bnd["low"] = 21.5
    bnd["high"] = 30.0
    bnd["bandpass"] = ""
    bnd["NET"] = 176.0
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 3.5
    bnd["A"] = 0.17752
    bnd["C"] = 0.76850
    bnd["NET_corr"] = 1.04
    bnd["pwv_poly"] = 0.973260, 0.082784, 0.001613
    bands["SAT_f030"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 38.75
    bnd["low"] = 30.0
    bnd["high"] = 47.5
    bnd["bandpass"] = ""
    bnd["NET"] = 217.6
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 3.5
    bnd["A"] = 0.41275
    bnd["C"] = 0.46189
    bnd["NET_corr"] = 1.01
    bnd["pwv_poly"] = 0.993708, 0.019467, 0.000420
    bands["SAT_f040"] = bnd
    
    bnd = OrderedDict()
    bnd["center"] = 25.75
    bnd["low"] = 21.5
    bnd["high"] = 30.0
    bnd["bandpass"] = ""
    bnd["NET"] = 176.0
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 3.5
    bnd["A"] = 0.17752
    bnd["C"] = 0.76850
    bnd["NET_corr"] = 1.04
    bnd["pwv_poly"] = 0.973260, 0.082784, 0.001613
    bands["CHSAT_f030"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 38.75
    bnd["low"] = 30.0
    bnd["high"] = 47.5
    bnd["bandpass"] = ""
    bnd["NET"] = 217.6
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 3.5
    bnd["A"] = 0.41275
    bnd["C"] = 0.46189
    bnd["NET_corr"] = 1.01
    bnd["pwv_poly"] = 0.993708, 0.019467, 0.000420
    bands["CHSAT_f040"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 91.5
    bnd["low"] = 77.0
    bnd["high"] = 106.0
    bnd["bandpass"] = ""
    bnd["NET"] = 278.1
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 3.5
    bnd["A"] = 0.19264
    bnd["C"] = 0.74873
    bnd["NET_corr"] = 1.15
    bnd["pwv_poly"] = 0.930223, 0.068510, 0.001771
    bands["CHLAT_f090"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 148.5
    bnd["low"] = 128.0
    bnd["high"] = 169.0
    bnd["bandpass"] = ""
    bnd["NET"] = 309.6
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 3.5
    bnd["A"] = 0.23072
    bnd["C"] = 0.69900
    bnd["NET_corr"] = 1.01
    bnd["pwv_poly"] = 0.763643, 0.231453, 0.006616
    bands["CHLAT_f150"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 91.5
    bnd["low"] = 77.0
    bnd["high"] = 106.0
    bnd["bandpass"] = ""
    bnd["NET"] = 287.0
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 3.5
    bnd["A"] = 0.27837
    bnd["C"] = 0.63703
    bnd["NET_corr"] = 1.15
    bnd["pwv_poly"] = 0.974286, 0.079438, 0.002080
    bands["SPLAT_f090"] = bnd
    
    bnd = OrderedDict()
    bnd["center"] = 148.5
    bnd["low"] = 128.0
    bnd["high"] = 169.0
    bnd["bandpass"] = ""
    bnd["NET"] = 268.8
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 3.5
    bnd["A"] = 0.22972
    bnd["C"] = 0.70032
    bnd["NET_corr"] = 1.01
    bnd["pwv_poly"] = 0.892955, 0.330207, 0.010174
    bands["SPLAT_f150"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 85.0
    bnd["low"] = 74.8
    bnd["high"] = 95.2
    bnd["bandpass"] = ""
    bnd["NET"] = 313.7
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 3.5
    bnd["A"] = 0.37419
    bnd["C"] = 0.51222
    bnd["NET_corr"] = 1.03
    bnd["pwv_poly"] = 0.980092, 0.061447, 0.001782
    bands["SAT_f085"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 145.0
    bnd["low"] = 129.1
    bnd["high"] = 161.0
    bnd["bandpass"] = ""
    bnd["NET"] = 335.3
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 3.5
    bnd["A"] = 0.19900
    bnd["C"] = 0.74041
    bnd["NET_corr"] = 1.00
    bnd["pwv_poly"] = 0.926762, 0.225619, 0.007903
    bands["SAT_f145"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 95.0
    bnd["low"] = 83.6
    bnd["high"] = 106.4
    bnd["bandpass"] = ""
    bnd["NET"] = 274.9
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 3.5
    bnd["A"] = 0.29310
    bnd["C"] = 0.61778
    bnd["NET_corr"] = 1.02
    bnd["pwv_poly"] = 0.972330, 0.085416, 0.002437
    bands["SAT_f095"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 155.0
    bnd["low"] = 138.0
    bnd["high"] = 172.1
    bnd["bandpass"] = ""
    bnd["NET"] = 359.7
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 3.5
    bnd["A"] = 0.18843
    bnd["C"] = 0.75417
    bnd["NET_corr"] = 1.00
    bnd["pwv_poly"] = 0.896329, 0.319420, 0.011034
    bands["SAT_f155"] = bnd
    
    bnd = OrderedDict()
    bnd["center"] = 85.0
    bnd["low"] = 74.8
    bnd["high"] = 95.2
    bnd["bandpass"] = ""
    bnd["NET"] = 313.7
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 3.5
    bnd["A"] = 0.37419
    bnd["C"] = 0.51222
    bnd["NET_corr"] = 1.03
    bnd["pwv_poly"] = 0.980092, 0.061447, 0.001782
    bands["CHSAT_f085"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 145.0
    bnd["low"] = 129.1
    bnd["high"] = 161.0
    bnd["bandpass"] = ""
    bnd["NET"] = 335.3
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 3.5
    bnd["A"] = 0.19900
    bnd["C"] = 0.74041
    bnd["NET_corr"] = 1.00
    bnd["pwv_poly"] = 0.926762, 0.225619, 0.007903
    bands["CHSAT_f145"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 95.0
    bnd["low"] = 83.6
    bnd["high"] = 106.4
    bnd["bandpass"] = ""
    bnd["NET"] = 274.9
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 3.5
    bnd["A"] = 0.29310
    bnd["C"] = 0.61778
    bnd["NET_corr"] = 1.02
    bnd["pwv_poly"] = 0.972330, 0.085416, 0.002437
    bands["CHSAT_f095"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 155.0
    bnd["low"] = 138.0
    bnd["high"] = 172.1
    bnd["bandpass"] = ""
    bnd["NET"] = 359.7
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 3.5
    bnd["A"] = 0.18843
    bnd["C"] = 0.75417
    bnd["NET_corr"] = 1.00
    bnd["pwv_poly"] = 0.896329, 0.319420, 0.011034
    bands["CHSAT_f155"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 227.0
    bnd["low"] = 197.9
    bnd["high"] = 256.1
    bnd["bandpass"] = ""
    bnd["NET"] = 683.2
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 3.5
    bnd["A"] = 0.33523
    bnd["C"] = 0.56279
    bnd["NET_corr"] = 1.01
    bnd["pwv_poly"] = 0.639310, 0.344960, 0.018401
    bands["CHLAT_f220"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 285.5
    bnd["low"] = 256.0
    bnd["high"] = 315.0
    bnd["bandpass"] = ""
    bnd["NET"] = 1679.7
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 3.5
    bnd["A"] = 0.40706
    bnd["C"] = 0.46930
    bnd["NET_corr"] = 1.01
    bnd["pwv_poly"] = 0.553007, 0.419329, 0.031032
    bands["CHLAT_f280"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 227.0
    bnd["low"] = 197.9
    bnd["high"] = 256.1
    bnd["bandpass"] = ""
    bnd["NET"] = 549.3
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 3.5
    bnd["A"] = 0.29109
    bnd["C"] = 0.62028
    bnd["NET_corr"] = 1.01
    bnd["pwv_poly"] = 0.825060, 0.534043, 0.034081
    bands["SPLAT_f220"] = bnd
    
    bnd = OrderedDict()
    bnd["center"] = 285.5
    bnd["low"] = 256.0
    bnd["high"] = 315.0
    bnd["bandpass"] = ""
    bnd["NET"] = 1295.9
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 3.5
    bnd["A"] = 0.33985
    bnd["C"] = 0.55678
    bnd["NET_corr"] = 1.01
    bnd["pwv_poly"] = 0.776403, 0.676221, 0.063374
    bands["SPLAT_f280"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 227.0
    bnd["low"] = 197.9
    bnd["high"] = 256.1
    bnd["bandpass"] = ""
    bnd["NET"] = 726.9
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 3.5
    bnd["A"] = 0.22747
    bnd["C"] = 0.70328
    bnd["NET_corr"] = 1.00
    bnd["pwv_poly"] = 0.867567, 0.403839, 0.027177
    bands["SAT_f220"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 285.5
    bnd["low"] = 256.0
    bnd["high"] = 315.0
    bnd["bandpass"] = ""
    bnd["NET"] = 1747.2
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 3.5
    bnd["A"] = 0.26587
    bnd["C"] = 0.65327
    bnd["NET_corr"] = 1.00
    bnd["pwv_poly"] = 0.830191, 0.512951, 0.050000
    bands["SAT_f280"] = bnd
    
    bnd = OrderedDict()
    bnd["center"] = 227.0
    bnd["low"] = 197.9
    bnd["high"] = 256.1
    bnd["bandpass"] = ""
    bnd["NET"] = 726.9
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 3.5
    bnd["A"] = 0.22747
    bnd["C"] = 0.70328
    bnd["NET_corr"] = 1.00
    bnd["pwv_poly"] = 0.867567, 0.403839, 0.027177
    bands["CHSAT_f220"] = bnd

    bnd = OrderedDict()
    bnd["center"] = 285.5
    bnd["low"] = 256.0
    bnd["high"] = 315.0
    bnd["bandpass"] = ""
    bnd["NET"] = 1747.2
    bnd["fknee"] = 50.0
    bnd["fmin"] = 0.01
    bnd["alpha"] = 3.5
    bnd["A"] = 0.26587
    bnd["C"] = 0.65327
    bnd["NET_corr"] = 1.00
    bnd["pwv_poly"] = 0.830191, 0.512951, 0.050000
    bands["CHSAT_f280"] = bnd

    cnf["bands"] = bands

    wafers = OrderedDict()

    wtypes = ["SPLAT_ULF", "CHLAT_LF", "SPLAT_LF", "SAT_LF", "CHSAT_LF", "CHLAT_MF", "SPLAT_MF", "SAT_MFL", "SAT_MFH", "CHSAT_MFL", "CHSAT_MFH", "CHLAT_HF", "SPLAT_HF", "SAT_HF", "CHSAT_HF"]
    # partial wafers will be counted as individual wafers since we make a full wafer then kill pixels
    wcnt = {
        "SPLAT_ULF": 4,
        "CHLAT_LF": 8*2,
        "SPLAT_LF": 9,
        "SAT_LF": (12) * 2,
        "CHSAT_LF": (12) * 2,
        "CHLAT_MF": 54*2,
        "SPLAT_MF": 54,
        "SAT_MFL": (12) * 6,
        "SAT_MFH": (12) * 6,
        "CHSAT_MFL": (12) * 6,
        "CHSAT_MFH": (12) * 6,
        "CHLAT_HF": 23*2,
        "SPLAT_HF": 18,
        "SAT_HF": (6 + 6) * 4,
        "CHSAT_HF": (6 + 6) * 4,
    }
    wpac = {
        "SPLAT_ULF": "RP",
        "CHLAT_LF": "RP",
        "SPLAT_LF": "RP",
        "SAT_LF": "RP",
        "CHSAT_LF": "RP",
        "CHLAT_MF": "RP",
        "SPLAT_MF": "RP",
        "SAT_MFL": "RP",
        "SAT_MFH": "HP",
        "CHSAT_MFL": "RP",
        "CHSAT_MFH": "HP",
        "CHLAT_HF": "HP",
        "SPLAT_HF": "HP",
        "SAT_HF": "HP",
        "CHSAT_HF": "HP",
    }
    wnp = {
        "SPLAT_ULF": 27,
        "CHLAT_LF": 48,
        "SPLAT_LF": 48,
        "SAT_LF": 12,
        "CHSAT_LF": 12,
        "CHLAT_MF": 432,
        "SPLAT_MF": 432,
        "SAT_MFL": 147,
        "SAT_MFH": 169,
        "CHSAT_MFL": 147,
        "CHSAT_MFH": 169,
        "CHLAT_HF": 469,
        "SPLAT_HF": 469,
        "SAT_HF": 469,
        "CHSAT_HF": 469,
    }
    wpixmm = {
        "SPLAT_ULF": 21.1,
        "CHLAT_LF": 16.1,
        "SPLAT_LF": 16.1,
        "SAT_LF": 31.1,
        "CHSAT_LF": 31.1,
        "CHLAT_MF": 5.3,
        "SPLAT_MF": 5.3,
        "SAT_MFL": 9.5,
        "SAT_MFH": 8.94,
        "CHSAT_MFL": 9.5,
        "CHSAT_MFH": 8.94,
        "CHLAT_HF": 5.2,
        "SPLAT_HF": 5.2,
        "SAT_HF": 5.2,
        "CHSAT_HF": 5.2,
    }
    wrhombgap = {
        "SPLAT_ULF": 2.827,
        "CHLAT_LF": 2.157,
        "SPLAT_LF": 2.157,
        "SAT_LF": 4.167,
        "CHSAT_LF": 4.167,
        "CHLAT_MF": 0.71,
        "SPLAT_MF": 0.71,
        "SAT_MFL": 1.273,
        "SAT_MFH": 0.71,
        "CHSAT_MFL": 1.273,
        "CHSAT_MFH": 0.71,
        "CHLAT_HF": 0.71,
        "SPLAT_HF": 0.71,
        "SAT_HF": 0.71,
        "CHSAT_HF": 0.71,
    }
    wbd = {
        "SPLAT_ULF": ["SPLAT_f020"],
        "CHLAT_LF": ["CHLAT_f030", "CHLAT_f040"],
        "SPLAT_LF": ["SPLAT_f030", "SPLAT_f040"],
        "SAT_LF": ["SAT_f030", "SAT_f040"],
        "CHSAT_LF": ["CHSAT_f030", "CHSAT_f040"],
        "CHLAT_MF": ["CHLAT_f090", "CHLAT_f150"],
        "SPLAT_MF": ["SPLAT_f090", "SPLAT_f150"],
        "SAT_MFL": ["SAT_f085", "SAT_f145"],
        "SAT_MFH": ["SAT_f095", "SAT_f155"],
        "CHSAT_MFL": ["CHSAT_f085", "CHSAT_f145"],
        "CHSAT_MFH": ["CHSAT_f095", "CHSAT_f155"],
        "CHLAT_HF": ["CHLAT_f220", "CHLAT_f280"],
        "SPLAT_HF": ["SPLAT_f220", "SPLAT_f280"],
        "SAT_HF": ["SAT_f220", "SAT_f280"],
        "CHSAT_HF": ["CHSAT_f220", "CHSAT_f280"],
    }
    windx = 0
    cardindx = 0
    for wt in wtypes:
        for ct in range(wcnt[wt]):
            wn = "{:02d}".format(windx)
            wf = OrderedDict()
            wf["type"] = wt
            wf["packing"] = wpac[wt]
            wf["rhombusgap"] = wrhombgap[wt]
            wf["npixel"] = wnp[wt]
            wf["pixsize"] = wpixmm[wt]
            wf["bands"] = wbd[wt]
            wf["card"] = "{:02d}".format(cardindx)
            cardindx += 1
            wafers[wn] = wf
            windx += 1

    cnf["wafers"] = wafers

    tubes = OrderedDict()

    woff = {
        "SPLAT_ULF": 0,
        "CHLAT_LF": 0,
        "SPLAT_LF": 0,
        "SAT_LF": 0,
        "CHSAT_LF": 0,
        "CHLAT_MF": 0,
        "SPLAT_MF": 0,
        "SAT_MFL": 0,
        "SAT_MFH": 0,
        "CHSAT_MFL": 0,
        "CHSAT_MFH": 0,
        "CHLAT_HF": 0,
        "SPLAT_HF": 0,
        "SAT_HF": 0,
        "CHSAT_HF": 0,
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
        "CHLAT_HF",
        "CHLAT_HF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
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
        "CHLAT_LF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_LF",
        "CHLAT_LF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_LF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_LF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_LF",
        "CHLAT_LF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_LF",
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
        "CHLAT_HF",
        "CHLAT_HF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
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
        "CHLAT_LF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_LF",
        "CHLAT_LF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_LF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_LF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_LF",
        "CHLAT_LF",
        "CHLAT_MF",
        "CHLAT_MF",
        "CHLAT_LF",
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
    ltubepos = [
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
        # tw is the wafer number in the tube. Here we use 6 for the 3 full and 3 partial wafers/tube
        for tw in range(1):
            off = 0
            for w, props in cnf["wafers"].items():
                if props["type"] == ttyp:
                    if off == woff[ttyp]:
                        tb["wafers"].append(w)
                        woff[ttyp] += 1
                        break
                    off += 1
        tb["platescale"] = 0.0047619
        tb["location"] = ltubepos[tindx]
        tubes[nm] = tb

    stubes = [
        "SAT_MFL",
        "SAT_MFH",
        "SAT_HF",
        "SAT_MFL",
        "SAT_MFH",
        "SAT_HF",
        "SAT_MFL",
        "SAT_MFH",
        "SAT_HF",
        "SAT_MFL",
        "SAT_MFH",
        "SAT_HF",
        "SAT_MFL",
        "SAT_MFH",
        "SAT_LF",
        "SAT_MFL",
        "SAT_MFH",
        "SAT_LF",
    ]
    stubepos = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    for tindx in range(18):
        nm = "ST{:d}".format(tindx)
        ttyp = stubes[tindx]
        tb = OrderedDict()
        tb["type"] = ttyp
        tb["waferspace"] = 124.
        tb["wafers"] = list()
        # HF tubes have 8 full wafers + 2 partial, all others 11+2
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
            # 29.0 deg/(422mm)
            tb["platescale"] = 0.0687
            tb["FOV_cut"] = 40.0
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
            # 29.4/(420mm)
            tb["platescale"] = 0.070
            tb["FOV_cut"] = 29.4
            tb["waferspace"] = 122.16
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
            # 29.4/(420mm)
            tb["platescale"] = 0.0701
            tb["FOV_cut"] = 29.4
            tb["waferspace"] = 121.85
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
            # 29.0/(490mm)
            tb["platescale"] = 0.0592
            tb["FOV_cut"] = 40.0
        tb["location"] = stubepos[tindx]
        tubes[nm] = tb

    #cnf["tubes"] = tubes
    
    chstubes = [
        "CHSAT_MFL",
        "CHSAT_MFH",
        "CHSAT_HF",
        "CHSAT_MFL",
        "CHSAT_MFH",
        "CHSAT_HF",
        "CHSAT_MFL",
        "CHSAT_MFH",
        "CHSAT_HF",
        "CHSAT_MFL",
        "CHSAT_MFH",
        "CHSAT_HF",
        "CHSAT_MFL",
        "CHSAT_MFH",
        "CHSAT_LF",
        "CHSAT_MFL",
        "CHSAT_MFH",
        "CHSAT_LF",
    ]
    chstubepos = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    for tindx in range(18):
        nm = "CHST{:d}".format(tindx)
        ttyp = chstubes[tindx]
        tb = OrderedDict()
        tb["type"] = ttyp
        tb["waferspace"] = 124.
        tb["wafers"] = list()
        # HF tubes have 8 full wafers + 2 partial, all others 11+2
        if ttyp == "CHSAT_HF":
            for tw in range(12):
                off = 0
                for w, props in cnf["wafers"].items():
                    if props["type"] == ttyp:
                        if off == woff[ttyp]:
                            tb["wafers"].append(w)
                            woff[ttyp] += 1
                            break
                        off += 1
            # 29.0 deg/(422mm)
            tb["platescale"] = 0.0687
            tb["FOV_cut"] = 40.0
        elif ttyp == "CHSAT_MFL":
            for tw in range(12):
                off = 0
                for w, props in cnf["wafers"].items():
                    if props["type"] == ttyp:
                        if off == woff[ttyp]:
                            tb["wafers"].append(w)
                            woff[ttyp] += 1
                            break
                        off += 1
            # 29.4/(420mm)
            tb["platescale"] = 0.070
            tb["FOV_cut"] = 29.4
            tb["waferspace"] = 122.16
        elif ttyp == "CHSAT_MFH":
            for tw in range(12):
                off = 0
                for w, props in cnf["wafers"].items():
                    if props["type"] == ttyp:
                        if off == woff[ttyp]:
                            tb["wafers"].append(w)
                            woff[ttyp] += 1
                            break
                        off += 1
            # 29.4/(420mm)
            tb["platescale"] = 0.0701
            tb["FOV_cut"] = 29.4
            tb["waferspace"] = 121.85
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
            # 29.0/(490mm)
            tb["platescale"] = 0.0592
            tb["FOV_cut"] = 40.0
        tb["location"] = chstubepos[tindx]
        tubes[nm] = tb

    cnf["tubes"] = tubes

    telescopes = OrderedDict()

    lfwhm = OrderedDict()
    lfwhm["SPLAT_f020"] = 11.4
    lfwhm["CHLAT_f030"] = 7.4
    lfwhm["SPLAT_f030"] = 8.4
    lfwhm["CHLAT_f040"] = 5.1
    lfwhm["SPLAT_f040"] = 5.8
    lfwhm["CHLAT_f090"] = 2.2
    lfwhm["SPLAT_f090"] = 2.5
    lfwhm["CHLAT_f150"] = 1.4
    lfwhm["SPLAT_f150"] = 1.6
    lfwhm["CHLAT_f220"] = 1.0
    lfwhm["SPLAT_f220"] = 1.1
    lfwhm["CHLAT_f280"] = 0.9
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
    tele["platescale"] = 0.0047619
    # This tube spacing in mm corresponds to 0.83 degrees projected on
    # the sky at a plate scale of 210 mm/deg or 0.0047619 deg/mm
    # The physical tube spacing is actually 210 mm.
    tele["tubespace"] =  174.3
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
    tele["platescale"] =  0.0047619
    # This tube spacing in mm corresponds to 0.83 degrees projected on
    # the sky at a plate scale of 210 mm/deg or 0.0047619 deg/mm
    # The physical tube spacing is actually 210 mm.
    tele["tubespace"] =  174.3
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
    tele["platescale"] = 0.0047619
    # This tube spacing in mm corresponds to 0.83 degrees projected on
    # the sky at a plate scale of 210 mm/deg or 0.0047619 deg/mm
    # The physical tube spacing is actually 210 mm.
    tele["tubespace"] =  174.3
    tele["fwhm"] = lfwhm
    telescopes["LAT2"] = tele

    #from the DSR
    sfwhm = OrderedDict()
    sfwhm["SAT_f030"] = 72.8
    sfwhm["SAT_f040"] = 72.8
    sfwhm["SAT_f085"] = 25.5
    sfwhm["SAT_f145"] = 25.5
    sfwhm["SAT_f095"] = 22.7
    sfwhm["SAT_f155"] = 22.7
    sfwhm["SAT_f220"] = 13.0
    sfwhm["SAT_f280"] = 13.0

    tele = OrderedDict()
    tele["tubes"] = ["ST0", "ST1", "ST2"]
    tele["platescale"] = 0.056689
    tele["tubespace"] = 700.0
    tele["fwhm"] = sfwhm
    telescopes["SAT0"] = tele

    tele = OrderedDict()
    tele["tubes"] = ["ST3", "ST4", "ST5"]
    tele["platescale"] = 0.056689
    tele["tubespace"] = 700.0
    tele["fwhm"] = sfwhm
    telescopes["SAT1"] = tele

    tele = OrderedDict()
    tele["tubes"] = ["ST6", "ST7", "ST8"]
    tele["platescale"] = 0.056689
    tele["tubespace"] = 700.0
    tele["fwhm"] = sfwhm
    telescopes["SAT2"] = tele

    tele = OrderedDict()
    tele["tubes"] = ["ST9", "ST10", "ST11"]
    tele["platescale"] = 0.056689
    tele["tubespace"] = 700.0
    tele["fwhm"] = sfwhm
    telescopes["SAT3"] = tele

    tele = OrderedDict()
    tele["tubes"] = ["ST12", "ST13", "ST14"]
    tele["platescale"] = 0.056689
    tele["tubespace"] = 700.0
    tele["fwhm"] = sfwhm
    telescopes["SAT4"] = tele

    tele = OrderedDict()
    tele["tubes"] = ["ST15", "ST16", "ST17"]
    tele["platescale"] = 0.056689
    tele["tubespace"] = 700.0
    tele["fwhm"] = sfwhm
    telescopes["SAT5"] = tele
    
    #from the DSR
    chsfwhm = OrderedDict()
    chsfwhm["CHSAT_f030"] = 72.8
    chsfwhm["CHSAT_f040"] = 72.8
    chsfwhm["CHSAT_f085"] = 25.5
    chsfwhm["CHSAT_f145"] = 25.5
    chsfwhm["CHSAT_f095"] = 22.7
    chsfwhm["CHSAT_f155"] = 22.7
    chsfwhm["CHSAT_f220"] = 13.0
    chsfwhm["CHSAT_f280"] = 13.0

    tele = OrderedDict()
    tele["tubes"] = ["CHST0", "CHST1", "CHST2"]
    tele["platescale"] = 0.056689
    tele["tubespace"] = 700.0
    tele["fwhm"] = chsfwhm
    telescopes["CHSAT0"] = tele

    tele = OrderedDict()
    tele["tubes"] = ["CHST3", "CHST4", "CHST5"]
    tele["platescale"] = 0.056689
    tele["tubespace"] = 700.0
    tele["fwhm"] = chsfwhm
    telescopes["CHSAT1"] = tele

    tele = OrderedDict()
    tele["tubes"] = ["CHST6", "CHST7", "CHST8"]
    tele["platescale"] = 0.056689
    tele["tubespace"] = 700.0
    tele["fwhm"] = chsfwhm
    telescopes["CHSAT2"] = tele

    tele = OrderedDict()
    tele["tubes"] = ["CHST9", "CHST10", "CHST11"]
    tele["platescale"] = 0.056689
    tele["tubespace"] = 700.0
    tele["fwhm"] = chsfwhm
    telescopes["CHSAT3"] = tele

    tele = OrderedDict()
    tele["tubes"] = ["CHST12", "CHST13", "CHST14"]
    tele["platescale"] = 0.056689
    tele["tubespace"] = 700.0
    tele["fwhm"] = chsfwhm
    telescopes["CHSAT4"] = tele

    tele = OrderedDict()
    tele["tubes"] = ["CHST15", "CHST16", "CHST17"]
    tele["platescale"] = 0.056689
    tele["tubespace"] = 700.0
    tele["fwhm"] = chsfwhm
    telescopes["CHSAT5"] = tele

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

    pl = ["A", "B"]
    hand = ["L", "R"]
    bandarr=["CHLAT_f030","CHLAT_f040"]

    dets = OrderedDict()
    for d in range(4):
        dprops = OrderedDict()
        dprops["wafer"] = "07"
        dprops["ID"] = d
        dprops["pixel"] = "000"
        bindx = d % 2
        dprops["band"] = bandarr[bindx]
        dprops["fwhm"] = 1.0
        dprops["pol"] = pl[bindx]
        dprops["handed"] = None
        dprops["card"] = "07"
        dprops["channel"] = d
        dprops["coax"] = 0
        dprops["bias"] = 0
        dprops["quat"] = np.array([0.0, 0.0, 0.0, 1.0])
        dname = "{}_{}_{}_{}".format("07", "000", dprops["band"], dprops["pol"])
        dets[dname] = dprops

    cnf["detectors"] = dets

    hw = Hardware()
    hw.data = cnf

    return hw
