"""
Functions providing the ensemble model and likelihood functions
"""
import numpy as np

# beast imports
# from beast.physicsmodel.prior_weights_dust import PriorWeightsDust

__all__ = ["lnlike", "lnprior", "lnprob"]


# temporary for development - remove and use PriorWeightsDust only
def _lognorm(x, max_pos, sigma=0.5, N=1.0):
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
    sqrt_2pi = 1.0 / np.sqrt(2 * np.pi)
    mu = np.log(max_pos) + sigma ** 2

    # avoid zero or negative due to log
    (indxs,) = np.where(x > 0)

    lnorm = np.zeros(len(x))

    log_x = np.log(x[indxs])
    normalization = sqrt_2pi / (x[indxs] * sigma)

    lnorm[indxs] = N * normalization * np.exp(-0.5 * ((log_x - mu) / sigma) ** 2)
    return lnorm


def _two_lognorm(xs, max_pos1, max_pos2, sigma1=0.5, sigma2=0.5, N1=1.0, N2=1.0):
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
    pointwise = _lognorm(xs, max_pos1, sigma=sigma1, N=N1) + _lognorm(
        xs, max_pos2, sigma=sigma2, N=N2
    )
    return pointwise


def lnlike(phi, beast_dust_priors, lnp_data, lnp_grid_vals):
    """
    Compute the log(likelihood) for the ensemble parameters

    Parameters
    ----------
    phi: floats
       ensemble parameters

    beast_dust_priors: PriorWeightsDust object
       contains the data and functions for the dust ensemble model

    lnp_data: dictonary
       contains arrays of the lnp values and indexs to the BEAST model grid

    lnp_grid_vals: dictonary
       contains arrays of the beast parameters and priors for the sparse
       lnp saved model grid points

    Returns
    -------
    log(likelihood): float
    """
    # unpack ensemble parameters (Av only)
    max_pos1, max_pos2, sigma1, sigma2, N12_ratio = phi

    # compute the ensemble model for all the model grid points for all stars
    #   temp code for development
    #   will change to using the PriorWeightsDust for production
    n_lnps, n_stars = beast_dust_priors.av_vals.shape
    new_prior = np.empty(beast_dust_priors.av_vals.shape, dtype=float)
    for k in range(n_stars):
        new_prior[:, k] = _two_lognorm(
            beast_dust_priors.av_vals[:, k],
            max_pos1,
            max_pos2,
            sigma1=sigma1,
            sigma2=sigma2,
            N1=1 - 1 / (N12_ratio + 1),
            N2=1 / (N12_ratio + 1),
        )
        if not np.isfinite(np.sum(new_prior[:, k])):
            print(new_prior[:, k])
            exit()

    # weights are those that adjust the saved likelihoods for the new
    # ensemble model (ensemble "priors")
    #   save as log to allow easy summing later
    weight_ratio = new_prior / beast_dust_priors.av_priors

    # compute the each star's integrated probability that it fits the new model
    # including the completeness function
    star_probs = np.sum(
        weight_ratio * lnp_grid_vals["completeness"] * np.exp(lnp_data["vals"]), axis=0
    )

    # remove any results that have zero integrated probabilities
    #   need to check why this is the case - all A(V) values of 0?
    (indxs,) = np.where(star_probs > 0.0)

    # print(np.sum(new_prior,axis=0))
    # print(np.sum(np.exp(lnp_data['vals']),axis=0))
    # print(np.sum(lnp_grid_vals['completeness'],axis=0))
    # print(star_probs)
    # print('---')
    # print(star_probs[indxs])
    # print('---')
    # print(np.sum(np.log(star_probs[indxs])))
    # print('===')
    # exit()

    # return the log product of the stars' probabilities
    return np.sum(np.log(star_probs[indxs]))


def lnprior(phi):
    """
    Compute the log(priors) for the ensemble parameters

    Parameters
    ----------
    phi: floats
       ensemble parameters

    Returns
    -------
    log(prior): floats
       0 if allowed
       -infinite if not allowed
    """
    # unpack ensemble parameters (Av only)
    max_pos1, max_pos2, sigma1, sigma2, N12_ratio = phi

    if (
        0.05 <= sigma1 < 2
        and 0.05 <= sigma2 < 2
        and 0 <= max_pos1 < 2
        and 0 <= max_pos2 < 3
        and max_pos1 < max_pos2
        and 0.01 < N12_ratio < 100
    ):
        return 0.0
    else:
        return -np.inf


def lnprob(phi, beast_dust_priors, lnp_data, lnp_grid_vals):
    """
    Compute the log(likelihood) for the ensemble parameters

    Parameters
    ----------
    phi: floats
       ensemble parameters

    beast_dust_priors: PriorWeightsDust object
       contains the data and functions for the dust ensemble model

    lnp_data: dictonary
       contains arrays of the lnp values and indexs to the BEAST model grid

    lnp_grid_vals: dictonary
       contains arrays of the beast parameters and priors for the sparse
       lnp saved model grid points

    Returns
    -------
    log(likelihood): float
    """
    ln_prior = lnprior(phi)

    if not np.isfinite(ln_prior):
        return -np.inf
    return ln_prior + lnlike(phi, beast_dust_priors, lnp_data, lnp_grid_vals)
