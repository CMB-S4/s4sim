# CMB-S4 Package s4sim

This package contains simulation tools that will be used for modeling of instrument
properties useful for trade studies and forecasts.


## Package Installation

This package uses setuptools for package installation.  You can use the setup.py
directly for installation or use pip.  Using pip makes it possible to uninstall the
package cleanly.  Installation to a specific location can be done with:

    %>  pip3 install . --prefix=/some/place

OR

    %>  python3 setup.py install --prefix /some/place

If you specify a `--prefix` option, then that location should already be in your PATH /
PYTHONPATH environment.  For example, you must have previously set up your environment
with:

    %>  export PATH="/some/place/bin:${PATH}"
    %>  export PYTHONPATH="/some/place/lib/python3.x/site-packages:${PYTHONPATH}"

Where `3.x` is the actual python version you are using.  If you do not specify a
`--prefix` option, the package will be installed to your python installation prefix.
**If** you are using a virtualenv or a conda environment for your development work
(always a good idea!) then it is safe to leave out the `--prefix` option and install the
package into your environment.

If you have used pip to do the installation, then you can uninstall with:

    %>  pip3 uninstall s4sim

If you want to be able to edit the source and have those changes show up in your
environment without re-installing, then you can use the package in develop / editable
mode.  This is done slightly differently depending on whether you use setup.py directly
or pip:

    %>  pip3 install -e . --prefix=/some/place

OR

    %>  python3 setup.py develop --prefix /some/place

These commands will fail unless the prefix is in your PYTHONPATH.


## Running Unit Tests

The unit tests are contained within the package itself.  This has the benefit that we
are certain to be running the tests on the installed package, rather than running tests
against a locally built (and perhaps different) version of the package.  After
installing the package (or using the develop / editable option), run tests with:

    %> python3 -c 'import s4sim.tests; s4sim.tests.run()'


## Package Development

When developing the code, you should install (e.g. with pip or conda) the "black" python
formatting package.  When making changes to the source, and before committing, run the
`format_source.sh` script in the top directory.  This applies formatting rules to all
the source files.

See DocDB XXXX for relevant coding guidelines.

## Hardware Configuration Package

The baseline hardware configuration is stored in s4sim/s4sim/hardware/ and is called via the get_example function.

• config.py defines the basic structure of the experiment including the detector and readout properties, the wafer properties, the tube properties, and the telescope properties

• sim.py defines the wafer layout and distribution of wafers in each optics tube

s4sim/s4sim/hardware/ also stored two outdated configurations (1shooter and 3shooter) that may be useful in future applications but are not the baseline design. The baseline design is called by the package init file, while these alternative configurations are not.

• 1shooter: This configuration has each SAT optics tube as its own telescope (i.e. all have independent scanning)

• 3shooter: This configuration has each SAT with three tubes (i.e. no independent scanning of optics tubes) but the tubes are not focused on the same patch of sky as they are in the baseline design. This could be adapted in the future if there are slight offsets on the area of the sky that different tubes see.
