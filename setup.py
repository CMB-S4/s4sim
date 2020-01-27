#!/usr/bin/env python3

import os
import sys
import re

import unittest
import shutil
import glob

from setuptools import find_packages, setup, Extension
from setuptools.command.test import test as TestCommand
from distutils.command.clean import clean as CleanCommand

import versioneer

# The setup options
setup_opts = dict()

# Entry points / scripts.  Add scripts here and define the main() of each
# script in s4sim.scripts.<foo>.main()
setup_opts["entry_points"] = {
    "console_scripts": ["s4sim_testscript = s4sim.scripts.testscript:main"]
}

setup_opts["name"] = "s4sim"
setup_opts["provides"] = ["s4sim"]
setup_opts["version"] = versioneer.get_version()
setup_opts["description"] = "CMB-S4 Simulation Package"
setup_opts["author"] = "CMB-S4 Collaboration"
setup_opts["author_email"] = "datamanagement@cmb-s4.org"
setup_opts["url"] = "https://github.com/CMB-S4/s4sim"
setup_opts["packages"] = find_packages(where=".")
setup_opts["license"] = "BSD"
setup_opts["requires"] = ["Python (>3.4.0)"]

# Command Class dictionary.
# Begin with the versioneer command class dictionary.
cmdcls = versioneer.get_cmdclass()

# Override the "clean" command to purge all build products.


class CustomClean(CleanCommand):
    def run(self):
        super().run()
        clean_files = [
            "./build",
            "./dist",
            "./__pycache__",
            "s4sim/__pycache__",
            "s4sim/tests/__pycache__",
            "s4sim/scripts/__pycache__",
            "./*.egg-info",
        ]
        for cf in clean_files:
            # Make paths absolute and relative to this path
            apaths = glob.glob(os.path.abspath(cf))
            for path in apaths:
                if os.path.isdir(path):
                    shutil.rmtree(path)
                elif os.path.isfile(path):
                    os.remove(path)
        return


cmdcls["clean"] = CustomClean

# Add command class to setup options
setup_opts["cmdclass"] = cmdcls

# Do the setup.
setup(**setup_opts)
