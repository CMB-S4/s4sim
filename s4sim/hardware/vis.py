# Copyright (c) 2020-2020 CMB-S4 Collaboration.
# Full license can be found in the top level "LICENSE" file.
"""Hardware visualization tools.
"""

import numpy as np

import toast
import toast.qarray as qa


default_band_colors = {
    "SPLAT_f020": (0.4, 0.4, 1.0, 0.2),
    "CHLAT_f030": (0.4, 0.4, 1.0, 0.2),
    "CHLAT_f040": (0.4, 0.4, 1.0, 0.2),
    "SPLAT_f030": (0.4, 0.4, 1.0, 0.2),
    "SPLAT_f040": (0.4, 0.4, 1.0, 0.2),
    "SAT_f030": (0.4, 0.4, 1.0, 0.2),
    "SAT_f040": (0.4, 0.4, 1.0, 0.2),
    "CHLAT_f090": (0.4, 1.0, 0.4, 0.2),
    "CHLAT_f150": (0.4, 1.0, 0.4, 0.2),
    "SPLAT_f090": (0.4, 1.0, 0.4, 0.2),
    "SPLAT_f150": (0.4, 1.0, 0.4, 0.2),
    "SAT_f085": (0.4, 1.0, 0.4, 0.2),
    "SAT_f145": (0.4, 1.0, 0.4, 0.2),
    "SAT_f095": (0.4, 1.0, 0.4, 0.2),
    "SAT_f155": (0.4, 1.0, 0.4, 0.2),
    "CHLAT_f220": (1.0, 0.4, 0.4, 0.2),
    "CHLAT_f280": (1.0, 0.4, 0.4, 0.2),
    "SPLAT_f220": (1.0, 0.4, 0.4, 0.2),
    "SPLAT_f280": (1.0, 0.4, 0.4, 0.2),
    "SAT_f220": (1.0, 0.4, 0.4, 0.2),
    "SAT_f280": (1.0, 0.4, 0.4, 0.2),
}


def plot_detectors(
    dets, outfile, width=None, height=None, labels=False, bandcolor=None, xieta=False,
):
    """Visualize a dictionary of detectors.

    This makes a simple plot of the detector positions on the projected
    focalplane.  The size of detector circles are controlled by the detector
    "fwhm" key, which is in arcminutes.  If the bandcolor is specified it will
    override the defaults.

    Args:
        outfile (str): Output PDF path.
        dets (dict): Dictionary of detector properties.
        width (float): Width of plot in degrees (None = autoscale).
        height (float): Height of plot in degrees (None = autoscale).
        labels (bool): If True, label each detector.
        bandcolor (dict, optional): Dictionary of color values for each band.
        xieta (bool):  If True, plot in Xi / Eta / Gamma coordinates rather
            than focalplane X / Y / Z.

    Returns:
        None

    """
    import matplotlib

    matplotlib.use("pdf")
    import matplotlib.pyplot as plt

    xaxis = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    zaxis = np.array([0.0, 0.0, 1.0], dtype=np.float64)

    n_det = len(dets)
    detnames = list(dets.keys())
    quats = np.array(
        [dets[detnames[x]]["quat"] for x in range(n_det)], dtype=np.float64
    )

    # arc_factor returns a scaling that can be used to reproject (X,
    # Y) to units corresponding to the angle subtended from (0,0,1) to
    # (X, Y, sqrt(X^2 + Y^2)).  This is called "ARC" (Zenithal
    # Equidistant) projection in FITS.
    def arc_factor(x, y):
        r = (x**2 + y**2) ** 0.5
        if r < 1e-6:
            return 1.0 + r**2 / 6
        return np.arcsin(r) / r

    if xieta:
        # Plotting as seen from observer in Xi / Eta / Gamma
        xi, eta, gamma = toast.instrument_coords.quat_to_xieta(quats)
        # xi, eta, gamma = quat.decompose_xieta(quats)

        detx = {
            detnames[k]: xi[k] * 180.0 / np.pi * arc_factor(xi[k], eta[k])
            for k in range(n_det)
        }
        dety = {
            detnames[k]: eta[k] * 180.0 / np.pi * arc_factor(xi[k], eta[k])
            for k in range(n_det)
        }

        for didx, dn in enumerate(detnames):
            wf = dets[dn]["wafer"]
            px = int(dets[dn]["pixel"])
            if px == 0:
                q = quats[didx]
                dd = qa.rotate(q, zaxis)

        # In Xi / Eta coordinates, gamma is measured clockwise from line of
        # decreasing elevation.  Here we convert into visualization X/Y
        # coordinatates measured counter clockwise from the X axis.
        polangs = {detnames[k]: 1.5 * np.pi - gamma[k] for k in range(n_det)}
    else:
        # Plotting in focalplane X / Y / Z coordinates.
        # Compute direction and orientation vectors
        dir = qa.rotate(quats, zaxis)
        orient = qa.rotate(quats, xaxis)

        small = np.fabs(1.0 - dir[:, 2]) < 1.0e-12
        not_small = np.logical_not(small)
        xp = np.zeros(n_det, dtype=np.float64)
        yp = np.zeros(n_det, dtype=np.float64)

        mag = np.arccos(dir[not_small, 2])
        ang = np.arctan2(dir[not_small, 1], dir[not_small, 0])
        xp[not_small] = mag * np.cos(ang)
        yp[not_small] = mag * np.sin(ang)

        polangs = {
            detnames[k]: np.arctan2(orient[k, 1], orient[k, 0]) for k in range(n_det)
        }

        detx = {
            detnames[k]: xp[k] * 180.0 / np.pi * arc_factor(xp[k], yp[k])
            for k in range(n_det)
        }
        dety = {
            detnames[k]: yp[k] * 180.0 / np.pi * arc_factor(xp[k], yp[k])
            for k in range(n_det)
        }

    if (width is None) or (height is None):
        # We are autoscaling.  Compute the angular extent of all detectors
        # and add some buffer.
        _y = np.array(list(dety.values()))
        _x = np.array(list(detx.values()))
        wmin, wmax = _x.min(), _x.max()
        hmin, hmax = _y.min(), _y.max()
        wbuf = 0.2 * (wmax - wmin)
        hbuf = 0.2 * (hmax - hmin)
        wmin -= wbuf
        wmax += wbuf
        hmin -= hbuf
        hmax += hbuf
        width = wmax - wmin
        height = hmax - hmin
        half_width = 0.5 * width
        half_height = 0.5 * height
    else:
        half_width = 0.5 * width
        half_height = 0.5 * height
        wmin = -half_width
        wmax = half_width
        hmin = -half_height
        hmax = half_height

    wafer_centers = dict()
    if labels:
        # We are plotting labels and will want to plot a wafer label for each
        # wafer.  To decide where to place the label, we find the average location
        # of all detectors from each wafer and put the label there.
        for d, props in dets.items():
            dwslot = props["wafer"]
            if dwslot not in wafer_centers:
                wafer_centers[dwslot] = []
            wafer_centers[dwslot].append((detx[d], dety[d]))
        for k in wafer_centers.keys():
            center = np.mean(wafer_centers[k], axis=0)
            size = (np.array(wafer_centers[k]) - center).std()
            wafer_centers[k] = (center, size)

    if bandcolor is None:
        bandcolor = default_band_colors
    wfigsize = 10.0
    hfigsize = wfigsize * (height / width)
    figdpi = 75
    hfigpix = int(figdpi * hfigsize)
    hpixperdeg = hfigpix / height

    fig = plt.figure(figsize=(wfigsize, hfigsize), dpi=figdpi)
    ax = fig.add_subplot(1, 1, 1)

    ax.set_xlim([wmin, wmax])
    ax.set_ylim([hmin, hmax])
    if xieta:
        ax.set_xlabel(r"Boresight $\xi$ Degrees", fontsize="large")
        ax.set_ylabel(r"Boresight $\eta$ Degrees", fontsize="large")
    else:
        ax.set_xlabel("Boresight X Degrees", fontsize="large")
        ax.set_ylabel("Boresight Y Degrees", fontsize="large")

    # Draw wafer labels in the background
    if labels:
        # The font size depends on the wafer size ... but keep it
        # between (0.01 and 0.1) times the size of the figure.
        for k, (center, size) in wafer_centers.items():
            fontpix = np.clip(0.7 * size * hpixperdeg, 0.01 * hfigpix, 0.10 * hfigpix)
            if fontpix < 1.0:
                fontpix = 1.0
            ax.text(
                center[0],
                center[1] + fontpix / hpixperdeg,
                k,
                color="k",
                fontsize=fontpix,
                horizontalalignment="center",
                verticalalignment="center",
                zorder=100,
                bbox=dict(fc="white", ec="none", pad=0.2, alpha=1.0),
            )

    for d, props in dets.items():
        band = props["band"]
        pixel = props["pixel"]
        pol = props["pol"]
        quat = np.array(props["quat"]).astype(np.float64)
        fwhm = props["fwhm"]

        # radius in degrees
        detradius = 0.5 * fwhm / 60.0

        # Position and polarization angle
        xpos, ypos = detx[d], dety[d]
        polang = polangs[d]

        detface = bandcolor[band]

        circ = plt.Circle(
            (xpos, ypos),
            radius=detradius,
            fc=detface,
            ec="black",
            linewidth=0.05 * detradius,
        )
        ax.add_artist(circ)

        ascale = 1.5

        xtail = xpos - ascale * detradius * np.cos(polang)
        ytail = ypos - ascale * detradius * np.sin(polang)
        dx = ascale * 2.0 * detradius * np.cos(polang)
        dy = ascale * 2.0 * detradius * np.sin(polang)

        detcolor = "black"
        if pol == "A":
            detcolor = (1.0, 0.0, 0.0, 1.0)
        if pol == "B":
            detcolor = (0.0, 0.0, 1.0, 1.0)

        ax.arrow(
            xtail,
            ytail,
            dx,
            dy,
            width=0.1 * detradius,
            head_width=0.3 * detradius,
            head_length=0.3 * detradius,
            fc=detcolor,
            ec="none",
            length_includes_head=True,
        )

        if labels:
            # Compute the font size to use for detector labels
            fontpix = 0.1 * detradius * hpixperdeg
            if fontpix < 1.0:
                fontpix = 1.0
            ax.text(
                xpos,
                ypos,
                pixel,
                color="k",
                fontsize=fontpix,
                horizontalalignment="center",
                verticalalignment="center",
                bbox=dict(fc="white", ec="none", pad=0.2, alpha=1.0),
            )
            labeloff = fontpix * len(pol) / hpixperdeg
            if dy < 0:
                labeloff = -labeloff
            ax.text(
                (xtail + 1.0 * dx + labeloff),
                (ytail + 1.0 * dy),
                pol,
                color="k",
                fontsize=fontpix,
                horizontalalignment="center",
                verticalalignment="center",
                bbox=dict(fc="none", ec="none", pad=0, alpha=1.0),
            )

    # Draw a "mini" coordinate axes for reference
    shortest = min(half_width, half_height)
    xmini = -0.7 * half_width
    ymini = -0.7 * half_height
    xlen = 0.06 * shortest
    ylen = 0.06 * shortest
    mini_width = 0.005 * shortest
    mini_head_width = 3 * mini_width
    mini_head_len = 3 * mini_width
    if xieta:
        aprops = [
            (xlen, 0, "-", r"$\xi$"),
            (0, ylen, "-", r"$\eta$"),
            (-xlen, 0, "--", "Y"),
            (0, -ylen, "--", "X"),
        ]
    else:
        aprops = [
            (xlen, 0, "-", "X"),
            (0, ylen, "-", "Y"),
            (-xlen, 0, "--", r"$\eta$"),
            (0, -ylen, "--", r"$\xi$"),
        ]
    for ap in aprops:
        lx = xmini + 1.5 * ap[0]
        ly = ymini + 1.5 * ap[1]
        lw = figdpi / 200.0
        ax.arrow(
            xmini,
            ymini,
            ap[0],
            ap[1],
            width=mini_width,
            head_width=mini_head_width,
            head_length=mini_head_len,
            fc="k",
            ec="k",
            linestyle=ap[2],
            linewidth=lw,
            length_includes_head=True,
        )
        ax.text(
            lx,
            ly,
            ap[3],
            color="k",
            fontsize=int(figdpi / 10),
            horizontalalignment="center",
            verticalalignment="center",
        )

    st = f"Focalplane Looking Towards Observer"
    if xieta:
        st = f"Focalplane on Sky From Observer"
    fig.suptitle(st)

    plt.savefig(outfile)
    plt.close()
    return








    # wmin = 1.0
    # wmax = -1.0
    # hmin = 1.0
    # hmax = -1.0
    # if (width is None) or (height is None):
    #     # We are autoscaling.  Compute the angular extent of all detectors
    #     # and add some buffer.
    #     for d, props in dets.items():
    #         quat = np.array(props["quat"]).astype(np.float64)
    #         dir = qa.rotate(quat, zaxis).flatten()
    #         if dir[0] > wmax:
    #             wmax = dir[0]
    #         if dir[0] < wmin:
    #             wmin = dir[0]
    #         if dir[1] > hmax:
    #             hmax = dir[1]
    #         if dir[1] < hmin:
    #             hmin = dir[1]
    #     wmin = np.arcsin(wmin) * 180.0 / np.pi
    #     wmax = np.arcsin(wmax) * 180.0 / np.pi
    #     hmin = np.arcsin(hmin) * 180.0 / np.pi
    #     hmax = np.arcsin(hmax) * 180.0 / np.pi
    #     wbuf = 0.1 * (wmax - wmin)
    #     hbuf = 0.1 * (hmax - hmin)
    #     wmin -= wbuf
    #     wmax += wbuf
    #     hmin -= hbuf
    #     hmax += hbuf
    #     width = wmax - wmin
    #     height = hmax - hmin
    # else:
    #     half_width = 0.5 * width
    #     half_height = 0.5 * height
    #     wmin = -half_width
    #     wmax = half_width
    #     hmin = -half_height
    #     hmax = half_height

    # if bandcolor is None:
    #     bandcolor = default_band_colors
    # xfigsize = 10.0
    # yfigsize = xfigsize * (height / width)
    # figdpi = 75
    # yfigpix = int(figdpi * yfigsize)
    # ypixperdeg = yfigpix / height

    # fig = plt.figure(figsize=(xfigsize, yfigsize), dpi=figdpi)
    # ax = fig.add_subplot(1, 1, 1)

    # ax.set_xlabel("Degrees", fontsize="large")
    # ax.set_ylabel("Degrees", fontsize="large")
    # ax.set_xlim([wmin, wmax])
    # ax.set_ylim([hmin, hmax])

    # for d, props in dets.items():
    #     band = props["band"]
    #     pixel = props["pixel"]
    #     pol = props["pol"]
    #     quat = np.array(props["quat"]).astype(np.float64)
    #     fwhm = props["fwhm"]

    #     # radius in degrees
    #     detradius = 0.5 * fwhm / 60.0

    #     # rotation from boresight
    #     rdir = qa.rotate(quat, zaxis).flatten()
    #     ang = np.arctan2(rdir[1], rdir[0])

    #     orient = qa.rotate(quat, xaxis).flatten()
    #     polang = np.arctan2(orient[1], orient[0])

    #     mag = np.arccos(rdir[2]) * 180.0 / np.pi
    #     xpos = mag * np.cos(ang)
    #     ypos = mag * np.sin(ang)

    #     detface = bandcolor[band]

    #     circ = plt.Circle(
    #         (xpos, ypos),
    #         radius=detradius,
    #         fc=detface,
    #         ec="black",
    #         linewidth=0.05 * detradius,
    #     )
    #     ax.add_artist(circ)

    #     ascale = 1.5

    #     xtail = xpos - ascale * detradius * np.cos(polang)
    #     ytail = ypos - ascale * detradius * np.sin(polang)
    #     dx = ascale * 2.0 * detradius * np.cos(polang)
    #     dy = ascale * 2.0 * detradius * np.sin(polang)

    #     detcolor = "black"
    #     if pol == "A":
    #         detcolor = (1.0, 0.0, 0.0, 1.0)
    #     if pol == "B":
    #         detcolor = (0.0, 0.0, 1.0, 1.0)

    #     ax.arrow(
    #         xtail,
    #         ytail,
    #         dx,
    #         dy,
    #         width=0.1 * detradius,
    #         head_width=0.3 * detradius,
    #         head_length=0.3 * detradius,
    #         fc=detcolor,
    #         ec="none",
    #         length_includes_head=True,
    #     )

    #     if labels:
    #         # Compute the font size to use for detector labels
    #         fontpix = 0.1 * detradius * ypixperdeg
    #         ax.text(
    #             (xpos),
    #             (ypos),
    #             pixel,
    #             color="k",
    #             fontsize=fontpix,
    #             horizontalalignment="center",
    #             verticalalignment="center",
    #             bbox=dict(fc="white", ec="none", pad=0.2, alpha=1.0),
    #         )
    #         xsgn = 1.0
    #         if dx < 0.0:
    #             xsgn = -1.0
    #         labeloff = 1.0 * xsgn * fontpix * len(pol) / ypixperdeg
    #         ax.text(
    #             (xtail + 1.0 * dx + labeloff),
    #             (ytail + 1.0 * dy),
    #             pol,
    #             color="k",
    #             fontsize=fontpix,
    #             horizontalalignment="center",
    #             verticalalignment="center",
    #             bbox=dict(fc="none", ec="none", pad=0, alpha=1.0),
    #         )

    # plt.savefig(outfile)
    # plt.close()
    # return


class clr:
    WHITE = "\033[97m"
    PURPLE = "\033[95m"
    BLUE = "\033[94m"
    GREEN = "\033[92m"
    YELLOW = "\033[93m"
    RED = "\033[91m"
    ENDC = "\033[0m"

    def disable(self):
        self.WHITE = ""
        self.PURPLE = ""
        self.BLUE = ""
        self.GREEN = ""
        self.YELLOW = ""
        self.RED = ""
        self.ENDC = ""


def summary_text(hw):
    """Print a textual summary of a hardware configuration.

    Args:
        hw (Hardware): A hardware dictionary.

    Returns:
        None

    """
    for obj, props in hw.data.items():
        nsub = len(props)
        print(
            "{}{:<12}: {}{:5d} objects{}".format(
                clr.WHITE, obj, clr.RED, nsub, clr.ENDC
            )
        )
        if nsub <= 2000:
            line = ""
            for k in list(props.keys()):
                if (len(line) + len(k)) > 72:
                    print("    {}{}{}".format(clr.BLUE, line, clr.ENDC))
                    line = ""
                line = "{}{}, ".format(line, k)
            if len(line) > 0:
                print("    {}{}{}".format(clr.BLUE, line.rstrip(", "), clr.ENDC))
        else:
            # Too many to print!
            print("    {}(Too many to print){}".format(clr.BLUE, clr.ENDC))

    return
