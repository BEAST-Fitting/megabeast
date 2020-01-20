#!/usr/bin/env python
"""
Script to run the MegaBEAST on the PHAT-like data.
"""

# system imports
from __future__ import absolute_import, division, print_function
import argparse

# MegaBEAST imports
from megabeast import megabeast
from megabeast import plot_input_data
from megabeast import parameter_maps


def megabeast_wrapper(
    megabeast_input_file, run_megabeast=False, diagnostic_plots=False, verbose=True
):
    """
    Wrapper to run the MegaBEAST.  This will eventually include making the diagnostic plots, too.

    Parameters
    ----------
    megabeast_input_file : string
        Name of the file that contains settings, filenames, etc

    run_megabeast : boolean (default=False)
        Run the MegaBEAST

    diagnostic_plots : boolean (default=False)
        Run code to create diagnostic plots (not implemented yet)

    verbose : boolean (default=True)
        print extra info

    """

    # run the MegaBEAST
    if run_megabeast:

        print("\n*********************")
        print("Running the MegaBEAST")
        print("*********************\n")

        megabeast.megabeast(megabeast_input_file, verbose=verbose)

        parameter_maps.parameter_maps(megabeast_input_file)

    # create diagnostic plots
    if diagnostic_plots:

        print("\n*************************")
        print("Creating diagnostic plots")
        print("*************************\n")

        plot_input_data.plot_input_data(megabeast_input_file, chi2_plot=[10])


if __name__ == "__main__":
    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-m", "--run_megabeast", help="Run the MegaBEAST", action="store_true"
    )
    parser.add_argument(
        "-p",
        "--diagnostic_plots",
        help="Generate diagnostic plots",
        action="store_true",
    )
    parser.add_argument("-v", "--verbose", help="Verbose output", action="store_true")
    parser.add_argument(
        "megabeast_input_file", help="path+filename for the MegaBEAST input file"
    )

    args = parser.parse_args()

    megabeast_wrapper(
        args.megabeast_input_file,
        run_megabeast=args.run_megabeast,
        diagnostic_plots=args.diagnostic_plots,
        verbose=args.verbose,
    )

    # print help if no arguments
    if not any(vars(args).values()):
        parser.print_help()
