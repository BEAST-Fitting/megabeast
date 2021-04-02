import numpy as np
import scipy.optimize as op

# from beast.physicsmodel.prior_weights_dust import PriorWeightsDust
from beast.physicsmodel.priormodel import PriorDustModel
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

        # format the physic model for the beast parameters
        #   uses the same format as the beast priors as these are equivalent to
        #   the megabeast physics model
        self.beast_params = ["Av", "Rv", "fA"]
        self.physics_model = {}
        for cparam in self.beast_params:
            cmod = self.fd_model[cparam]
            self.physics_model[cparam] = {"name": cmod["name"]}
            for cname, cval in zip(cmod["varnames"], cmod["varinit"]):
                self.physics_model[cparam][cname] = cval

    def start_params(self):
        """
        Get the start parameters for the fit

        Returns
        -------
        values, names : tuple
            names give the parameters names and values for all the submodels
            with non-fixed priors
        """
        dmod = self.fd_model
        names = []
        values = []
        for ckey in dmod.keys():
            if dmod[ckey]["prior"]["name"] != "fixed":
                for cparam, cval in zip(dmod[ckey]["varnames"], dmod[ckey]["varinit"]):
                    names.append(f"{ckey}_{cparam}")
                    values.append(cval)
        return (values, values)

    def lnlike(self, phi, beast_dust_priors, lnp_data, lnp_grid_vals):
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
        # update the values in the physics model using the input ensemble parameters
        k = 0
        for cparam in self.beast_params:
            cmod = self.fd_model[cparam]
            if cmod["prior"]["name"] != "fixed":
                for cvar in cmod["varnames"]:
                    self.physics_model[cparam][cvar] = phi[k]
                    k += 1
        # print(self.physics_model)

        # compute the ensemble model for all the model grid points for all stars
        n_lnps, n_stars = beast_dust_priors.av_vals.shape
        mb_physics_dust = np.zeros((n_lnps, n_stars), dtype=float)

        av_prior = PriorDustModel(self.physics_model["Av"])
        rv_prior = PriorDustModel(self.physics_model["Rv"])
        fa_prior = PriorDustModel(self.physics_model["fA"])

        gvals = np.isfinite(lnp_data["vals"])
        mb_physics_dust[gvals] = (
            av_prior(beast_dust_priors.av_vals[gvals])
            * rv_prior(beast_dust_priors.rv_vals[gvals])
            * fa_prior(beast_dust_priors.fA_vals[gvals])
        )

        # test for a spoiler stars
        if np.any(~np.isfinite(np.sum(mb_physics_dust, axis=0))):
            for k, psum in enumerate(np.sum(mb_physics_dust, axis=0)):
                if not np.isfinite(psum):
                    print(beast_dust_priors.av_vals[gvals, k])
                    print(mb_physics_dust[:, k])
                    raise ValueError(
                        f"star #{k} has a non-finite integrated prob for this model"
                    )

        # weights are those that adjust the saved likelihoods for the new
        # ensemble model (ensemble "priors")
        #   save as log to allow easy summing later
        weight_ratio = mb_physics_dust / beast_dust_priors.av_priors

        # compute the each star's integrated probability that it fits the new model
        # including the completeness function
        star_probs = np.sum(
            weight_ratio * lnp_grid_vals["completeness"] * np.exp(lnp_data["vals"]),
            axis=0,
        )

        for k in range(n_stars):
            if not np.isfinite(star_probs[k]):
                print(k)
                print(weight_ratio[:, k])
                print(lnp_grid_vals["completeness"][:, k])
                # print(lnp_data["vals"][:, k])
                exit()

        # remove any results that have zero integrated probabilities
        #   need to check why this is the case - all A(V) values of 0?
        (indxs,) = np.where(star_probs > 0.0)

        # print(" ")
        # print(star_probs)
        # print(np.sum(np.log(star_probs[indxs])))
        # print(" ")
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
        # check each dust parameter
        dmod = self.fd_model
        k = 0
        for cdust in dmod.keys():
            cprior = dmod[cdust]["prior"]
            cname = cprior["name"]
            if cname != "fixed":
                if cname == "flat":
                    for i, cparam in enumerate(dmod[cdust]["varnames"]):
                        if (
                            phi[k] < cprior["minmax"][i][0]
                            or phi[k] > cprior["minmax"][i][1]
                        ):
                            return -np.inf
                        k += 1
                else:
                    raise ValueError(f"{cname} prior not supported")
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

    if not np.isfinite(ln_prior):
        return -np.inf
    return ln_prior + megabeast_model.lnlike(
        phi, beast_dust_priors, lnp_data, lnp_grid_vals
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
    lnp_data = read_lnp_data(lnp_filename)

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
        megabeast_model.start_params()[0],
        args=(megabeast_model, beast_dust_priors, lnp_data, lnp_grid_vals),
        method="Nelder-Mead",
    )

    # next step would be to
    # run through MCMC to fully sample likelihood
    # maybe include option not to run MCMC

    # print("output")
    # print(megabeast_model.physics_model)
    print(result)
    return result["x"]
