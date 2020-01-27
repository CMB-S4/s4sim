# Copyright (c) 2020-2020 CMB-S4 Collaboration.
# Full license can be found in the top level "LICENSE" file.
"""Example script.
"""

import argparse

# Use relative imports for things within the package.
from ..fake import fake_function


def main():
    parser = argparse.ArgumentParser(
        description="This is an example script.",
        usage="s4sim_testscript [file, [file]] ...",
    )

    parser.add_argument("files", type=str, nargs="+", help="Input file(s)")

    args = parser.parse_args()

    for file in args.files:
        print("Input file {}...".format(file), flush=True)

    fake_function()

    return
