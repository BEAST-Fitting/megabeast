"""
Script to run the MegaBEAST on BEAST results.
"""

# system
import argparse
import os

# other packages
from tqdm import trange
import numpy as np
import scipy.optimize as op
from astropy.io import fits

# beast
from beast.physicsmodel.prior_weights_dust import PriorWeightsDust
from beast.tools import read_beast_data

# megabeast
from megabeast.read_megabeast_input import read_megabeast_input
from megabeast.ensemble_model import lnprob


def megabeast_single():
    """
    Run the MegaBEAST on a single set of BEAST results
    """
    pass


def megabeast_image(megabeast_input_file, verbose=True):
    """
    Run the MegaBEAST on an image of BEAST results.  The BEAST results
    are given as spatially-reordered BEAST outputs with a file of lnp results
    for each pixel in the image.

    Parameters
    ----------
    megabeast_input_file : string
        Name of the file that contains settings, filenames, etc

    verbose : boolean (default=True)
        print extra info

    """
    # read in the settings from the file
    mb_settings = read_megabeast_input(megabeast_input_file)

    # setup the megabeast model including defining the priors
    #   - dust distribution model
    #   - stellar populations model (later)

    # use nstars image to setup for each pixel
    nstars_image, nstars_header = fits.getdata(
        mb_settings["nstars_filename"], header=True
    )
    n_x, n_y = nstars_image.shape

    # read in the beast data that is needed by all the pixels
    beast_data = {}
    # - SED data
    beast_data.update(
        read_beast_data.read_sed_data(
            mb_settings["beast_seds_filename"], param_list=["Av"]  # , "Rv", "f_A"]
        )
    )
    # - max completeness
    beast_data.update(
        read_beast_data.read_noise_data(
            mb_settings["beast_noise_filename"], param_list=["completeness"],
        )
    )
    beast_data["completeness"] = np.max(beast_data["completeness"], axis=1)

    # setup for output
    pixel_fit_status = np.full((n_x, n_y), False, dtype=bool)
    n_fit_params = len(mb_settings["fit_param_names"])
    best_fit_images = np.zeros((n_x, n_y, n_fit_params), dtype=float) + np.nan

    # loop over the pixels with non-zero entries in the nstars image
    for i in trange(n_x, desc="x pixels"):
        for j in trange(n_y, desc="y pixels", leave=False):
            # for i in [6]:
            #    for j in [6]:
            if verbose:
                print("working on (%i,%i)" % (i, j))
            if nstars_image[i, j] >= mb_settings["min_for_fit"]:
                pixel_fit_status[i, j] = True
                # get the saved sparse likelihoods
                lnp_filename = mb_settings[
                    "lnp_file_prefix"
                ] + "_{0}_{1}_lnp.hd5".format(j, i)
                lnp_data = read_beast_data.read_lnp_data(
                    lnp_filename, nstars=nstars_image[i, j], shift_lnp=True,
                )

                # get the completeness and BEAST model parameters for the
                #   same grid points as the sparse likelihoods
                lnp_grid_vals = read_beast_data.get_lnp_grid_vals(beast_data, lnp_data)

                # initialize the ensemble model with the parameters used
                # for the saved BEAST model run results
                #   currently only dust parameters allowed
                #   for testing -> only Av
                avs = lnp_grid_vals["Av"]
                rvs = [3.1]  # beast_data['Rv']
                fAs = [1.0]  # beast_data['f_A']
                beast_dust_priors = PriorWeightsDust(
                    avs,
                    mb_settings["av_prior_model"],
                    rvs,
                    mb_settings["rv_prior_model"],
                    fAs,
                    mb_settings["fA_prior_model"],
                )

                # standard minimization to find initial values
                def chi2(args):
                    return -1.0 * lnprob(*args)
                result = op.minimize(
                    chi2,
                    [0.25, 2.0, 0.5, 0.5, 1],
                    args=(beast_dust_priors, lnp_data, lnp_grid_vals),
                    method="Nelder-Mead",
                )
                best_fit_images[i, j, :] = result["x"]
                # print(result)
                # print(result['x'])
                # print(result['success'])

                # then run through MCMC to fully sample likelihood
                #    include option not to run MCMC

    # output results
    #    - best fit
    #    - megabeast parameter 1D pPDFs
    #    - MCMC chain

    master_header = nstars_header
    # Now, write the maps to disk

    # check that the directory exists
    if not os.path.exists("./" + mb_settings["projectname"] + "_megabeast/"):
        os.makedirs("./" + mb_settings["projectname"] + "_megabeast/")

    for k, cname in enumerate(mb_settings["fit_param_names"]):

        hdu = fits.PrimaryHDU(best_fit_images[:, :, k], header=master_header)

        # Save to FITS file
        hdu.writeto(
            "%s_megabeast/%s_%s_bestfit.fits"
            % (mb_settings["projectname"], mb_settings["projectname"], cname),
            overwrite=True,
        )


if __name__ == "__main__":
    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "megabeast_input_file",
        help="Name of the file that contains settings, filenames, etc",
    )
    parser.add_argument("-v", "--verbose", help="verbose output", action="store_true")
    args = parser.parse_args()

    megabeast_image(args.megabeast_input_file, verbose=args.verbose)
