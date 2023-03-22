# Copyright (c) 2020-2020 CMB-S4 Collaboration.
# Full license can be found in the top level "LICENSE" file.

import pickle

import numpy as np

from toast.pipeline_tools import Telescope, Focalplane, Site, Schedule, CES
from toast.timing import function_timer, Timer
from toast.utils import Logger
import toast.qarray as qa

from .. import hardware


XAXIS, YAXIS, ZAXIS = np.eye(3)

# Note, this is for atmospheric sims only and doesn't affect the wafer/tube scaling to the sky.
# For SATs, this is FOV/2*(1+2/sqrt(3)), FOV UHF=35 deg, FOV others=29 deg
FOCALPLANE_RADII_DEG = {
    "LAT0": 3.9,
    "LAT1": 3.9,
    "LAT2": 3.9,
    "SAT0": 17.5,
    "SAT1": 17.5,
    "SAT2": 17.5,
    "SAT3": 17.5,
    "SAT4": 14.5,
    "SAT5": 14.5,
}


class S4Telescope(Telescope):
    def __init__(self, name, site=None):
        if (
            (site is None and (name == "LAT0"  or name == "LAT1"))
            or (site is not None and site.upper() == "CHILE")
        ):
            site = Site("Atacama", lat="-22.958064", lon="-67.786222", alt=5200)
            offset = 0
        elif site is None or site.upper() == "POLE":
            site = Site("Pole", lat="-89.991067", lon="-44.650000", alt=2843)
            offset = 1024
            # HACK begin
            if name == "LAT2":
                offset = 0
            # HACK end
        else:
            raise RuntimeError("Unknown site: {}".format(site))
        super().__init__(name, site=site)
        self.id = {
            # Use the same telescope index for telescopes of the same type
            # in the same place to re-use the atmospheric simulations
            #'LAT0' : 0, 'LAT1' : 1, 'LAT2' : 2, 'SAT0' : 3, 'SAT1' : 4...
            "LAT0": 1 + offset,
            "LAT1": 1 + offset,
            "LAT2": 1 + offset,
            "SAT0": 8 + offset,
            "SAT1": 8 + offset,
            "SAT2": 8 + offset,
            "SAT3": 8 + offset,
            "SAT4": 8 + offset,
            "SAT5": 8 + offset,
        }[name]


def add_hw_args(parser):
    parser.add_argument(
        "--hardware", required=False, default=None, help="Input hardware file"
    )
    parser.add_argument(
        "--site", required=False, default=None, help="Observing site: Chile or Pole"
    )
    parser.add_argument(
        "--thinfp",
        required=False,
        type=np.int,
        help="Thin the focalplane by this factor",
    )
    parser.add_argument(
        "--radiusfp-deg",
        required=False,
        type=np.float,
        help="Reject detectors outside of this radius from the boresight. "
        "To reject inside the radius, provide a negative value.",
    )
    parser.add_argument(
        "--bands",
        required=True,
        help="Comma-separated list of bands: SPLAT_f020 (20 GHz, SPLAT), "
        "CHLAT_f030 (30 GHz CHLAT), CHLAT_f040 (40 GHz, CHLAT), "
        "SPLAT_f030 (30 GHz SPLAT), SPLAT_f040 (40 GHz, SPLAT), "
        "SAT_f030 (30 GHz, SAT), SAT_f040 (40 GHz, SAT), CHLAT_f090 (90 GHz, CHLAT), "
        "CHLAT_f150 (150 GHz, CHLAT), SPLAT_f090 (90 GHz, SPLAT), "
        "SPLAT_f150 (150 GHz, SPLAT), SAT_f085 (85 GHz, SAT), "
        "SAT_f0145 (145 GHz, SAT), SAT_f095 (95 GHz, SAT), SAT_f155 (155 GHz, SAT), "
        "CHLAT_f220 (220 GHz, CHLAT), CHLAT_f280 (280 GHz, CHLAT), "
        "SPLAT_f220 (220 GHz, SPLAT), SPLAT_f280 (280 GHz, SPLAT), "
        "SAT_f220 (220 GHz, SAT), SAT_f280 (280 GHz, SAT)."
    )
    parser.add_argument(
        "--tubes",
        required=False,
        help="Comma-separated list of optics tubes:\n"
        "LAT0-CHLAT_LF : LT63, LT66, LT67, LT70, LT75, LT78, LT79, LT82. "
        "LAT0-CHLAT_MF : LT19..LT22, LT25..LT31, LT34..LT62, LT64, LT65, LT68, LT69, "
        "                LT71..LT74, LT76, LT77, LT80, LT81, LT83, LT84. "
        "LAT0-CHLAT_HF : LT0..LT18, LT23, LT24, LT32, LT33. "
        "LAT1-CHLAT_LF : LT148, LT151, LT152, LT155, LT160, LT163, LT164, LT167. "
        "LAT1-CHLAT_MF : LT104..LT107, LT110..LT116, LT119..LT150, LT153, LT154, "
        "                LT156..LT159, LT161, LT162, LT165, LT166, LT168, LT169. "
        "LAT1-CHLAT_HF : LT85..LT103, LT108, LT109, LT117, LT118. "
        "LAT2-SPLAT_ULF : LT178, LT182, LT184, LT188. "
        "LAT2-SPLAT_LF : LT232, LT234, LT236, LT239, LT242, LT245, LT248, LT251, LT254. "
        "LAT2-SPLAT_MF : LT170..LT176, LT180, LT186, LT189..LT206, LT208, LT210, LT212, LT214, LT216, LT218, LT220, LT222, LT224, LT226, LT228, LT230, LT231, LT233, LT235, LT237, LT238, LT240, LT241, LT243, LT244, LT246, LT247, LT249, LT250, LT252, LT253. "
        "LAT2-SPLAT_HF : LT177, LT179, LT181, LT183, LT185, LT187, LT207, LT209, LT211, LT213, LT215, LT217, LT219, LT221, LT223, LT225, LT227, LT229. "
        "SAT0-SAT_MFL : ST0. "
        "SAT0-SAT_MFH : ST1. "
        "SAT0-SAT_HF : ST2. "
        "SAT1-SAT_MFL : ST3. "
        "SAT1-SAT_MFH : ST4. "
        "SAT1-SAT_HF : ST5. "
        "SAT2-SAT_MFL : ST6. "
        "SAT2-SAT_MFH  : ST7. "
        "SAT2-SAT_HF : ST8. "
        "SAT3-SAT_MFL : ST9. "
        "SAT3-SAT_MFH  : ST10. "
        "SAT3-SAT_HF : ST11. "
        "SAT4-SAT_LF : ST14. "
        "SAT4-SAT_MFL : ST12. "
        "SAT4-SAT_MFH : ST13. "
        "SAT5-SAT_LF : ST17. "
        "SAT5-SAT_MFL : ST15. "
        "SAT6-SAT_MFH : ST16."
    )
    parser.add_argument(
        "--telescope",
        required=False,
        help="Telescope, one of: LAT0, LAT1, LAT2, SAT0, SAT1, SAT2, "
        "SAT3, SAT4, SAT5",
    )
    return


class BandParams:
    def __init__(self, band_name, band_data):
        self.name = band_name
        self.net = band_data["NET"] * 1e-6  # uK -> K
        self.fknee = band_data["fknee"] * 1e-3  # mHz -> Hz
        self.fmin = band_data["fmin"] * 1e-3  # mHz -> Hz
        self.alpha = band_data[band]["alpha"]
        # overwrite the hardware model value of 3.5, to a more realistic value
        #self.alpha = 1.0
        self.A = band_data["A"]
        self.C = band_data["C"]
        self.lower = band_data["low"]  # GHz
        self.center = band_data["center"]  # GHz
        self.upper = band_data["high"]  # GHz
        return


class DetectorParams:
    def __init__(self, det_data, band, wafer, tube, telescope, index):
        """
        Args:
            det_data (dict) :  Dictionary of detector parameters
                from s4sim.hardware
            band (BandParams) :  band parameters act as defaults
            wafer (int) :  wafer number
            tube (str) :  tube name
            telescope (str) :  telescope name
            index (int) :  RNG index
        """

        def get_par(key, default, scale=1):
            if key in det_data:
                return det_data[key] * scale
            else:
                return default

        self.band = band.name
        self.det_data = det_data
        self.net = get_par("NET", band.net, 1e-6)  # uK -> K
        self.fknee = get_par("fknee", band.fknee, 1e-3)  # mHz -> Hz
        self.fmin = get_par("fmin", band.fmin, 1e-3)  # mHz -> Hz
        self.alpha = get_par("alpha", band.alpha)
        # overwrite the hardware model value of 3.5, to a more realistic value
        self.alpha = 1
        self.A = get_par("A", band.A)
        self.C = get_par("C", band.C)
        self.lower = get_par("low", band.lower)  # GHz
        self.center = get_par("center", band.center)  # GHz
        self.upper = get_par("high", band.upper)  # GHz
        # ensure that the center frequency of band is center of upper and lower bands
        self.center = 0.5 * (self.lower + self.upper)
        self.width = self.upper - self.lower
        self.wafer = wafer
        self.tube = tube
        self.telescope = telescope
        self.index = index
        return

    def get_dict(self):
        det_dict = {
            "NET": self.net,
            "fknee": self.fknee,
            "fmin": self.fmin,
            "alpha": self.alpha,
            "A": self.A,
            "C": self.C,
            "quat": self.det_data["quat"],
            "fwhm": self.det_data["fwhm"],
            "freq": self.center,
            "bandcenter_ghz": self.center,
            "bandwidth_ghz": self.width,
            "index": self.index,
            "telescope": self.telescope,
            "tube": self.tube,
            "wafer": self.wafer,
            "band": self.band,
        }
        return det_dict


@function_timer
def get_hardware(args, comm, verbose=False):
    """ Get the hardware configuration, either from file or by simulating.
    Then trim it down to the bands that were selected.
    """
    log = Logger.get()
    timer = Timer()
    timer.start()
    telescope = get_telescope(args, comm, verbose=verbose)
    if comm.world_rank == 0:
        if args.hardware:
            log.info("Loading hardware configuration from {}...".format(args.hardware))
            if args.hardware.endswith(".pkl"):
                with open(args.hardware, "rb") as fin:
                    hw = pickle.load(fin)
            else:
                hw = hardware.Hardware(args.hardware)
            timer.report_clear("Load {}".format(args.hardware))
        else:
            log.info("Simulating default hardware configuration")
            hw = hardware.get_example()
            timer.report_clear("Get example hardware")
            hw.data["detectors"] = hardware.sim_telescope_detectors(hw, telescope.name)
            timer.report_clear("Get telescope detectors")
        # Construct a running index for all detectors across all
        # telescopes for independent noise realizations
        det_index = {}
        for idet, det in enumerate(sorted(hw.data["detectors"])):
            det_index[det] = idet
        match = {"band": args.bands.replace(",", "|")}
        if args.tubes is None:
            tubes = None
        else:
            tubes = args.tubes.split(",")
        # If one provides both telescopes and tubes, the tubes matching *either*
        # will be concatenated
        # hw = hw.select(telescopes=[telescope.name], tubes=tubes, match=match)
        hw = hw.select(tubes=tubes, match=match)
        ndetector = len(hw.data["detectors"])
        if ndetector == 0:
            raise RuntimeError(
                "No detectors match query: telescope={}, "
                "tubes={}, match={}".format(telescope.name, tubes, match)
            )
        if args.thinfp:
            # Only accept a fraction of the detectors for
            # testing and development
            thin_index = {}
            for idet, det in enumerate(sorted(hw.data["detectors"])):
                thin_index[det] = idet
            delete_detectors = []
            for det_name in hw.data["detectors"].keys():
                if (thin_index[det_name] // 2) % args.thinfp != 0:
                    delete_detectors.append(det_name)
            for det_name in delete_detectors:
                del hw.data["detectors"][det_name]
        if args.radiusfp_deg:
            # Only accept detectors inside the given radius
            radius = np.radians(args.radiusfp_deg)
            # Comparing to the cosine of the radius avoids repeated
            # arccos evaluations
            cosradius = np.cos(radius)
            inside = radius > 0
            delete_detectors = []
            for det_name, det_data in hw.data["detectors"].items():
                det_quat = det_data["quat"]
                vec = qa.rotate(det_quat, ZAXIS)
                rcos = np.dot(vec, ZAXIS)
                if (rcos < cosradius and inside) or (rcos > cosradius and not inside):
                    delete_detectors.append(det_name)
            for det_name in delete_detectors:
                del hw.data["detectors"][det_name]

        ndetector = len(hw.data["detectors"])
        log.info(
            "Telescope = {} tubes = {} bands = {}, thinfp = {} radiusfp = {} "
            "matches {} detectors".format(
                telescope.name,
                args.tubes,
                args.bands,
                args.thinfp,
                args.radiusfp_deg,
                ndetector,
            )
        )
        timer.report_clear("Trim detectors")
    else:
        hw = None
        det_index = None
    if comm.comm_world is not None:
        hw = comm.comm_world.bcast(hw)
        det_index = comm.comm_world.bcast(det_index)
    return hw, telescope, det_index


@function_timer
def get_telescope(args, comm, verbose=False):
    """ Determine which telescope matches the detector selections
    """
    if args.telescope is not None:
        return S4Telescope(args.telescope, site=args.site)
    telescope = None
    if comm.world_rank == 0:
        hwexample = hardware.get_example()
        tubes = args.tubes.split(",")
        for tube in tubes:
            for telescope_name, telescope_data in hwexample.data["telescopes"].items():
                if tube in telescope_data["tubes"]:
                    if telescope is None:
                        telescope = S4Telescope(telescope_name, site=args.site)
                    elif telescope.name != telescope_name:
                        raise RuntimeError(
                            "Tubes '{}' span more than one telescope".format(tubes)
                        )
                    break
            else:
                raise RuntimeError(
                    "Failed to match tube = '{}' with a telescope".format(tube)
                )
    if comm.comm_world is not None:
        telescope = comm.comm_world.bcast(telescope)
    return telescope


def get_focalplane(args, comm, hw, det_index, verbose=False):
    """ Translate hardware configuration into a TOAST focalplane dictionary
    """
    if comm.world_rank == 0:
        detector_data = {}
        band_params = {}
        for band_name, band_data in hw.data["bands"].items():
            band_params[band_name] = BandParams(band_name, band_data)
        # User may force the effective focal plane radius to be larger
        # than the default.  This will widen the size of the simulated
        # atmosphere but has no other effect for the time being.
        fpradius = None
        try:
            fpradius = args.focalplane_radius_deg
        except:
            pass
        if fpradius is None:
            fpradius = 0
        for det_name, det_data in hw.data["detectors"].items():
            # RNG index for this detector
            index = det_index[det_name]
            wafer = det_data["wafer"]
            # Determine which tube has this wafer
            for tube_name, tube_data in hw.data["tubes"].items():
                if wafer in tube_data["wafers"]:
                    break
            # Determine which telescope has this tube
            for telescope_name, telescope_data in hw.data["telescopes"].items():
                if tube_name in telescope_data["tubes"]:
                    break
            else:
                raise RuntimeError(
                    "Unable to match tube {} to a telescope".format(tube_name)
                )
            fpradius = max(fpradius, FOCALPLANE_RADII_DEG[telescope_name])
            det_params = DetectorParams(
                det_data,
                band_params[det_data["band"]],
                wafer,
                tube_name,
                telescope_name,
                index,
            )
            detector_data[det_name] = det_params.get_dict()
        # Create a focal plane in the telescope
        focalplane = Focalplane(
            detector_data=detector_data,
            sample_rate=args.sample_rate,
            radius_deg=fpradius,
        )
    else:
        focalplane = None
    if comm.comm_world is not None:
        focalplane = comm.comm_world.bcast(focalplane)
    return focalplane


@function_timer
def load_focalplanes(args, comm, schedules, verbose=False):
    """ Attach a focalplane to each of the schedules.

    Args:
        schedules (list) :  List of Schedule instances.
            Each schedule has two members, telescope
            and ceslist, a list of CES objects.
    Returns:
        detweights (dict) : Inverse variance noise weights for every
            detector across all focal planes. In [K_CMB^-2].
            They can be used to bin the TOD.
    """
    # log = Logger.get()
    timer = Timer()
    timer.start()

    # Load focalplane information

    timer1 = Timer()
    timer1.start()
    hw, telescope, det_index = get_hardware(args, comm, verbose=verbose)
    focalplane = get_focalplane(args, comm, hw, det_index, verbose=verbose)
    telescope.focalplane = focalplane

    if comm.world_rank == 0 and verbose:
        timer1.report_clear("Collect focaplane information")

    for schedule in schedules:
        # Replace the telescope created from reading the observing schedule but
        # keep the weather object
        weather = schedule.telescope.site.weather
        schedule.telescope = telescope
        schedule.telescope.site.weather = weather

    detweights = telescope.focalplane.detweights

    timer.stop()
    if (comm.comm_world is None or comm.world_rank == 0) and verbose:
        timer.report("Loading focalplane")
    return detweights
