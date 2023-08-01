# Copyright (c) 2020-2023 CMB-S4 Collaboration.
# Full license can be found in the top level "LICENSE" file.
"""Hardware models for use in analysis.

This module contains code for simulating a hardware model and
dumping / loading hardware information to / from disk.

"""

# These are simply namespace imports for convenience.

from .config import Hardware, sim_nominal

from .sim import (
    sim_detectors_toast,
    sim_detectors_physical_optics,
    sim_telescope_detectors,
)

from .vis import plot_detectors, summary_text
