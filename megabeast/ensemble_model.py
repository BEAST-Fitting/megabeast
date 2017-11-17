"""
Functions providing the ensemble model and likelihood functions
"""

# system imports
from __future__ import (absolute_import, division, print_function)

# other package imports
import numpy as np

# beast imports
from beast.physicsmodel.prior_weights_dust import PriorWeightsDust

# temporary for development - remove and use PriorWeightsDust only
def _lognorm(x, max_pos, sigma=0.5, N=1.):
    """
    Lognormal distribution

    Parameters
    ----------
    xs: vector
       x values

    max_pos: float
       Position of the lognormal function's maximum

    sigma: float
       Sigma of the lognormal function

    N: floats
       Multiplicative factor

    Returns
    -------
    lognormal computed on the x grid
    """
    sqrt_2pi = 1. / np.sqrt(2 * np.pi)
    mu = np.log(max_pos) + sigma**2

    # avoid zero or negative due to log
    indxs, = np.where(x > 0)

    lnorm = np.zeros(len(x))

    log_x = np.log(x[indxs])
    normalization = sqrt_2pi / (x[indxs] * sigma)
    
    lnorm[indxs] = (N * normalization 
                    * np.exp(-0.5 * ((log_x - mu) / sigma)**2))
    return lnorm

def _two_lognorm(xs, 
                 max_pos1, max_pos2, 
                 sigma1=0.5, sigma2=0.5, 
                 N1=1., N2=1.):
    """
    Mixture of 2 lognormal functions

    Parameters
    ----------
    xs: vector
       x values

    max_pos1: float
       Position of the lognormal function's maximum for component 1
    max_pos2: float
       Position of the lognormal function's maximum for component 2

    sigma1: float
       Sigma of the lognormal function for component 1
    sigma2: float
       Sigma of the lognormal function for component 2

    N1: floats
       Multiplicative factor for component 1
    N2: floats
       Multiplicative factor for component 2

    Returns
    -------
    Mixture model: (LOGNORM1 + LOGNORM2)
    """
    pointwise = (_lognorm(xs, max_pos1, sigma=sigma1, N=N1)
                 + _lognorm(xs, max_pos2, sigma=sigma2, N=N2))
    return pointwise

def lnlike(phi, model_weights, lnp_data, beast_on_lnp):
    """
    Compute the log(likelihood) for the ensemble parameters

    Parameters
    ----------
    phi: floats
       ensemble parameters

    model_weights: PriorWeightsDust object
       contains the data and functions for the dust ensemble model

    lnp_data: dictonary
       contains arrays of the lnp values and indexs to the BEAST model grid

    beast_on_lnp: dictonary
       contains arrays of the beast parameters and priors for the sparse
       lnp saved model grid points

    Returns
    -------
    log(likelihood): float
    """
    # unpack ensemble parameters (Av only)
    max_pos1, max_pos2, sigma1, sigma2, N1, N2 = phi

    # compute the ensemble model for all the model grid points for all the stars
    #   temp code for development
    #   will change to using the PriorWeightsDust for production
    n_lnps, n_stars = model_weights.av_vals.shape
    new_prior = np.empty(model_weights.av_vals.shape, dtype=float)
    for k in range(n_stars):
        new_prior[:,k] = _two_lognorm(model_weights.av_vals[:,k],
                                      max_pos1, max_pos2,
                                      sigma1=sigma1, sigma2=sigma2,
                                      N1=N1, N2=N2)
        if not np.isfinite(np.sum(new_prior[:,k])):
            print(new_prior[:,k])
            exit()
        #if np.sum(new_prior[:,k]) == 0.:
        #    print(new_prior[:,k])
        #    print(model_weights.av_vals[:,k])
        #    exit()
            

    # weights are those that adjust the saved likelihoods for the new
    # ensemble model (ensemble "priors")
    #   save as log to allow easy summing later
    weight_ratio = new_prior/beast_on_lnp['prior_weight']

    # compute the each star's integrated probability that it fits the new model
    # including the completeness function
    star_probs = np.sum(weight_ratio
                        *beast_on_lnp['completeness']
                        *np.exp(lnp_data['vals']),axis=0)

    # remove any results that have zero integrated probabilities
    #   need to check why this is the case - all A(V) values of 0?
    indxs, = np.where(star_probs > 0.)

    #print(np.sum(new_prior,axis=0))
    #print(np.sum(np.exp(lnp_data['vals']),axis=0))
    #print(np.sum(beast_on_lnp['completeness'],axis=0))
    #print(star_probs)
    #print('---')
    #print(star_probs[indxs])
    #print('---')
    #print(np.sum(np.log(star_probs[indxs])))
    #print('===')
    #exit()

    # return the log product of the stars' probabilities 
    return np.sum(np.log(star_probs[indxs]))

def lnprior(phi):
    """
    Compute the log(priors) for the ensemble parameters

    Parameters
    ----------
    phi: floats
       ensemble parameters

    Returns:
    --------
    log(prior): floats
       0         if allowed
       -infinite if not allowed
    """
    # unpack ensemble parameters (Av only)
    max_pos1, max_pos2, sigma1, sigma2, N1, N2 = phi

    if (0.05 <= sigma1 < 2
        and 0.05 <= sigma2 < 2
        and 0 <= max_pos1 < 2
        and 0 <= max_pos2 < 3
        and max_pos1 < max_pos2
        and N1 >= 0
        and N2 >= 0):
        return 0.0
    else:
        return -np.inf

def lnprob(phi, model_weights, lnp_data, beast_on_lnp):
    """
    Compute the log(likelihood) for the ensemble parameters

    Parameters
    ----------
    phi: floats
       ensemble parameters

    model_weights: PriorWeightsDust object
       contains the data and functions for the dust ensemble model

    lnp_data: dictonary
       contains arrays of the lnp values and indexs to the BEAST model grid

    beast_on_lnp: dictonary
       contains arrays of the beast parameters and priors for the sparse
       lnp saved model grid points

    Returns
    -------
    log(likelihood): float
    """
    ln_prior = lnprior(phi)

    if not np.isfinite(ln_prior):
        return -np.inf
    return ln_prior + lnlike(phi, model_weights, lnp_data, beast_on_lnp)
    
