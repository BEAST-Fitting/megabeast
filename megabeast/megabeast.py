"""
Script to run the MegaBEAST on BEAST results.
"""

# system
from __future__ import (absolute_import, division, print_function)
import argparse
import os

# other packages
from tqdm import (tqdm, trange)
import numpy as np
import scipy.optimize as op
from astropy.io import fits

# beast
from beast.physicsmodel.prior_weights_dust import PriorWeightsDust

# megabeast
from .read_megabeast_input import read_megabeast_input
from .beast_data import (read_beast_data, extract_beast_data,
                        read_lnp_data)
from .ensemble_model import lnprob

import pdb

def megabeast(megabeast_input_file, verbose=True):
    """
    Run the MegaBEAST on each of the spatially-reordered BEAST outputs.

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
    nstars_image, nstars_header = fits.getdata(mb_settings['nstars_filename'], header=True)
    n_x, n_y = nstars_image.shape

    # read in the beast data that is needed by all the pixels
    beast_data = read_beast_data(mb_settings['beast_seds_filename'],
                                 mb_settings['beast_noise_filename'],
                                 beast_params=['completeness',
                                               'Av'])#,'Rv','f_A'])

    # setup for output
    pixel_fit_status = np.full((n_x, n_y), False, dtype=bool)
    n_fit_params = len(mb_settings['fit_param_names'])
    best_fit_images = np.zeros((n_x, n_y, n_fit_params), dtype=float) + np.nan

    # loop over the pixels with non-zero entries in the nstars image
    for i in trange(n_x, desc="x pixels"):
        for j in trange(n_y, desc="y pixels", leave=False):
    #for i in [6]:
    #    for j in [6]:
            if verbose:
                print("working on (%i,%i)"%(i,j))
            if nstars_image[i,j] >= mb_settings['min_for_fit']:
                pixel_fit_status[i,j] = True
                # get the saved sparse likelihoods
                lnp_filename = mb_settings['lnp_file_prefix']+"_%i_%i_lnp.hd5"%(j, i)
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
                dustpriors = PriorWeightsDust(avs, mb_settings['av_prior_model'],
                                              rvs, mb_settings['rv_prior_model'],
                                              fAs, mb_settings['fA_prior_model'])

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

    # check that the directory exists
    if not os.path.exists('./'+mb_settings['projectname'] + '_megabeast/'):
        os.makedirs('./'+mb_settings['projectname'] + '_megabeast/')

    for k, cname in enumerate(mb_settings['fit_param_names']):

        hdu = fits.PrimaryHDU(best_fit_images[:,:,k], header=master_header)

        # Save to FITS file
        hdu.writeto("%s_megabeast/%s_%s_bestfit.fits"%(mb_settings['projectname'], mb_settings['projectname'], cname),
                    overwrite=True)

        

if __name__ == '__main__':
    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("megabeast_input_file",
                        help="Name of the file that contains settings, filenames, etc")
    parser.add_argument("-v", "--verbose", help="verbose output",
                        action="store_true")
    args = parser.parse_args()

    megabeast(args.megabeast_input_file, verbose=args.verbose)

