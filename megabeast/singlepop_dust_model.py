import numpy as np
import scipy.optimize as op

from beast.physicsmodel.prior_weights_dust import PriorWeightsDust
from beast.tools.read_beast_data import (
    read_lnp_data,
    get_lnp_grid_vals,
)

__all__ = ["MB_Model"]


class MB_Model:
    """
    MegaBEAST model that provides member functions to compute
    the likelihood and priors for a specific physical model
    """

    def __init__(self, params):
        self.fd_model = params.fd_model

    def start_params(self):
        """
        Return the start parameters for the fit
        """
        dmod = self.fd_model
        return np.array([dmod["Av_model"]["m_init"], dmod["Av_model"]["s_init"]])

    def lnlike(phi, megabeast_model, beast_dust_priors, lnp_data, lnp_grid_vals):
        """
        Compute the log(likelihood) for the ensemble parameters

        Parameters
        ----------
        phi: floats
           ensemble parameters

        megabeast_model : dict
            MegaBEAST physical model including priors

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
        mval, sval, N12_ratio = phi

        # compute the ensemble model for all the model grid points for all stars
        #   temp code for development
        #   will change to using the PriorWeightsDust for production
        n_lnps, n_stars = beast_dust_priors.av_vals.shape
        new_prior = np.zeros(beast_dust_priors.av_vals.shape, dtype=float)
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

    def lnprior(self, phi):
        """
        Compute the log(priors) for the ensemble parameters

        Parameters
        ----------
        phi: floats
           ensemble parameters

        megabeast_model : dict
            MegaBEAST physical model including priors

        Returns
        -------
        log(prior): floats
           0 if allowed
           -infinite if not allowed
        """
        # unpack ensemble parameters (Av only)
        mval, sval = phi
        print(mval, sval)

        # check each dust parameter
        dmod = self.fd_model
        for cdust in dmod.keys():
            cprior = dmod[cdust]["prior"]
            if cprior["name"] != "fixed":
                if cprior["name"] == "flat":
                    if mval < cprior["m_min"] or mval > cprior["m_max"]:
                        return -np.inf
                    if mval < cprior["s_min"] or mval > cprior["s_max"]:
                        return -np.inf
                else:  # prior not supported
                    exit()
                    # raise warning instead
        return 0.0


def lnprob(phi, megabeast_model, beast_dust_priors, lnp_data, lnp_grid_vals):
    """
    Compute the log(likelihood) for the ensemble parameters

    Parameters
    ----------
    phi: floats
       ensemble parameters

    megabeast_model : class
        MegaBEAST physical model including priors

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
    ln_prior = megabeast_model.lnprior(phi)
    print(ln_prior)
    exit()

    if not np.isfinite(ln_prior):
        return -np.inf
    return ln_prior + megabeast_model.lnlike(
        phi, megabeast_model, beast_dust_priors, lnp_data, lnp_grid_vals
    )


def fit_ensemble(beast_data, lnp_filename, beast_priormodel, megabeast_model):
    """
    Run the MegaBEAST on a single set of BEAST results.

    Parameters
    ----------
    beast_data : dict
        information about the BEAST runs including SED grid and noise model

    lnp_filename : string
        file with posteriors from BEAST fitting

    beast_priormodel : dict
        BEAST prior model information

    megabeast_model : dict
        MegaBEAST physical model including priors

    Returns
    -------
    fit_results : array
        set of best fit parameters
    """
    # get the saved sparse likelihoods
    lnp_data = read_lnp_data(lnp_filename, shift_lnp=True)

    # get the completeness and BEAST model parameters for the
    #   same grid points as the sparse likelihoods
    lnp_grid_vals = get_lnp_grid_vals(beast_data, lnp_data)

    # compute the BEAST prior weights
    #  needed so the BEAST posteriors updated with the MegaBEAST model
    # ***currently only AV ensemble model supported***
    avs = lnp_grid_vals["Av"]
    rvs = lnp_grid_vals["Rv"]
    fAs = lnp_grid_vals["f_A"]
    beast_dust_priors = PriorWeightsDust(
        avs,
        beast_priormodel["AV"],
        rvs,
        beast_priormodel["RV"],
        fAs,
        beast_priormodel["fA"],
    )

    # standard minimization to find initial values
    def chi2(*args):
        return -1.0 * lnprob(*args)

    result = op.minimize(
        chi2,
        megabeast_model.start_params(),
        args=(megabeast_model, beast_dust_priors, lnp_data, lnp_grid_vals),
        # method="Nelder-Mead",
    )

    # next step would be to
    # run through MCMC to fully sample likelihood
    # maybe include option not to run MCMC

    return result["x"]
