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
from beast.tools.read_beast_data import (
    read_lnp_data,
    get_lnp_grid_vals,
    read_sed_data,
    read_noise_data,
)

# megabeast
from megabeast.read_megabeast_input import read_megabeast_input
from megabeast.ensemble_model import lnprob


def fit_ensemble(beast_data, lnp_filename, beast_priormodel, nstars_expected=None):
    """
    Run the MegaBEAST on a single set of BEAST results.

    Parameters
    ----------
    beast_data : dict
        information about the BEAST runs including SED grid and noise model

    lnp_filename : string
        file with posteriors from BEAST fitting

    beast_priormodel : dict
        dictionary of the BEAST prior model information

    nstars_expected : int
        number of stars expected, used as a check

    Returns
    -------
    fit_results : array
        set of best fit parameters
    """
    # get the saved sparse likelihoods
    lnp_data = read_lnp_data(lnp_filename, nstars=nstars_expected, shift_lnp=True)

    # get the completeness and BEAST model parameters for the
    #   same grid points as the sparse likelihoods
    lnp_grid_vals = get_lnp_grid_vals(beast_data, lnp_data)

    # compute the BEAST prior weights
    #  needed so the BEAST posteriors updated with the MegaBEAST model
    # ***currently only AV ensemble model supported***
    avs = lnp_grid_vals["Av"]
    rvs = [3.1]  # beast_data['Rv']
    fAs = [1.0]  # beast_data['f_A']
    beast_dust_priors = PriorWeightsDust(
        avs,
        beast_priormodel["AV"],
        rvs,
        beast_priormodel["RV"],
        fAs,
        beast_priormodel["fA"],
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

    # next step would be to
    # run through MCMC to fully sample likelihood
    # maybe include option not to run MCMC

    return result["x"]


def megabeast_single():
    """
    Run the MegaBEAST on a single set of BEAST results.
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
    params = read_megabeast_input(megabeast_input_file)

    # use nstars image to setup for each pixel
    nstars_image, nstars_header = fits.getdata(params["nstars_filename"], header=True)
    n_x, n_y = nstars_image.shape

    # read in the beast data that is needed by all the pixels
    beast_data = {}
    # - SED data
    beast_data.update(read_sed_data(params["beast_seds_filename"], param_list=["Av"]))
    # - max completeness
    beast_data.update(
        read_noise_data(params["beast_noise_filename"], param_list=["completeness"],)
    )
    beast_data["completeness"] = np.max(beast_data["completeness"], axis=1)

    # BEAST prior model
    beast_pmodel = {}
    beast_pmodel["AV"] = params["av_prior_model"]
    beast_pmodel["RV"] = params["rv_prior_model"]
    beast_pmodel["fA"] = params["fA_prior_model"]

    # setup for output
    pixel_fit_status = np.full((n_x, n_y), False, dtype=bool)
    n_fit_params = len(params["fit_param_names"])
    best_fit_images = np.zeros((n_x, n_y, n_fit_params), dtype=float) + np.nan

    # loop over the pixels with non-zero entries in the nstars image
    for i in trange(n_x, desc="x pixels"):
        for j in trange(n_y, desc="y pixels", leave=False):

            if nstars_image[i, j] >= params["min_for_fit"]:
                pixel_fit_status[i, j] = True

                # filename with saved BEAST posteriors
                lnp_prefix = params["lnp_file_prefix"]
                lnp_filename = f"{lnp_prefix}_{j}_{i}_lnp.hd5"

                best_fit_params = fit_ensemble(
                    beast_data,
                    lnp_filename,
                    beast_pmodel,
                    nstars_expected=nstars_image[i, j],
                )

                best_fit_images[i, j, :] = best_fit_params

    # output results (* = future)
    #    - best fit
    #    - *megabeast parameter 1D pPDFs
    #    - *MCMC chain

    # Write the maps to disk
    master_header = nstars_header

    # check that the directory exists
    dpath = "./%s_megabeast/" % (params["projectname"])
    if not os.path.exists(dpath):
        os.makedirs(dpath)

    for k, cname in enumerate(params["fit_param_names"]):
        hdu = fits.PrimaryHDU(best_fit_images[:, :, k], header=master_header)
        hdu.writeto(
            "%s_megabeast/%s_%s_bestfit.fits"
            % (params["projectname"], params["projectname"], cname),
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
