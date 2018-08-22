"""
Script to run the MegaBEAST on BEAST results.
"""

# system
from __future__ import (absolute_import, division, print_function)
import argparse

# other packages
from tqdm import (tqdm, trange)
import numpy as np
import scipy.optimize as op
from astropy.io import fits

# beast
from beast.physicsmodel.prior_weights_dust import PriorWeightsDust

# megabeast
from beast_data import (read_beast_data, extract_beast_data,
                        read_lnp_data)
from ensemble_model import lnprob


def megabeast(projectname, min_for_fit=20, verbose=True):
    """
    Run the MegaBEAST on each of the spatially-reordered BEAST outputs.

    Parameters
    ----------
    projectname : string
        Project name to use (basename for files)

    min_for_fit : int (default=20)
        minimum number of stars need in a pixel for fit

    verbose : boolean (default=True)
        print extra info

    """

    # setup the megabeast model including defining the priors
    #   - dust distribution model
    #   - stellar populations model (later)

    # use nstars image to setup for each pixel
    nstars_filename = "%s_nstars.fits" % (projectname)
    nstars_image, nstars_header = fits.getdata(nstars_filename, header=True)
    n_x, n_y = nstars_image.shape

    # read in the beast data that is needed by all the pixels
    beast_seds_filename = "%s_seds.grid.hd5"%(projectname)
    beast_noise_filename = "%s_noisemodel.hd5"%(projectname)
    beast_data = read_beast_data(beast_seds_filename,
                                 beast_noise_filename,
                                 beast_params=['completeness',
                                               'Av'])#,'Rv','f_A'])

    # setup for output
    pixel_fit_status = np.full((n_x, n_y), False, dtype=bool)
    fit_param_names = ['Av1', 'Av2', 'sigma1', 'sigma2', 'N1', 'N2']
    n_fit_params = len(fit_param_names)
    best_fit_images = np.zeros((n_x, n_y, n_fit_params), dtype=float)

    # loop over the pixels with non-zero entries in the nstars image
    for i in trange(n_x, desc="x pixels"):
        for j in trange(n_y, desc="y pixels", leave=False):
    #for i in [6]:
    #    for j in [6]:
            if verbose:
                print("working on (%i,%i)"%(i,j))
            if nstars_image[i,j] >= min_for_fit:
                pixel_fit_status[i,j] = True
                # get the saved sparse likelihoods
                lnp_filename = "%s_%i_%i_lnp.hd5"%(projectname, j, i)
                lnp_data = read_lnp_data(lnp_filename, nstars_image[i,j])

                # get the completeness and BEAST model parameters for the
                #   same grid points as the sparse likelihoods
                beast_on_lnp = extract_beast_data(beast_data, lnp_data)

                # initialize the ensemble model with the parameters used
                # for the saved BEAST model run results
                #   currently only dust parameters allowed
                #   for testing -> only Av
                avs = beast_on_lnp['Av']
                rvs = [3.1]#beast_data['Rv']
                fAs = [1.0]#beast_data['f_A']
                av_prior_model = {'name': 'flat'}
                rv_prior_model = {'name': 'flat'}
                fA_prior_model = {'name': 'flat'}
                dustpriors = PriorWeightsDust(avs, av_prior_model,
                                              rvs, rv_prior_model,
                                              fAs, fA_prior_model)

                # standard minimization to find initial values
                chi2 = lambda * args: -1.0*lnprob(*args)
                N12_init = 0.5*nstars_image[i,j]
                result = op.minimize(chi2,
                                     [0.1,0.7,0.5,0.5,N12_init,N12_init],
                                     args=(dustpriors, lnp_data, beast_on_lnp),
                                     method='Nelder-Mead')
                best_fit_images[i,j,:] = result['x']
                #print(result)
                #print(result['x'])
                #print(result['success'])

                # then run through MCMC to fully sample likelihood
                #    include option not to run MCMC

    # output results
    #    - best fit
    #    - megabeast parameter 1D pPDFs
    #    - MCMC chain

    master_header = nstars_header
    # Now, write the maps to disk

    for k, cname in enumerate(fit_param_names):

        hdu = fits.PrimaryHDU(best_fit_images[:,:,k], header=master_header)

        # Save to FITS file
        hdu.writeto("%s_%s_bestfit.fits"%(projectname, cname),
                    overwrite=True)


if __name__ == '__main__':
    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("projectname",
                        help="project name to use (basename for files)")
    parser.add_argument("--min_for_fit", default=20, type=int,
                        help="minimum number of stars need in a pixel for fit")
    parser.add_argument("-v", "--verbose", help="verbose output",
                        action="store_true")
    args = parser.parse_args()

    megabeast(args.projectname, min_for_fit=args.min_for_fit, verbose=args.verbose)

