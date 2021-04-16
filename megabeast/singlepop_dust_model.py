import numpy as np
import scipy.optimize as op

# from beast.physicsmodel.prior_weights_dust import PriorWeightsDust
from beast.physicsmodel.priormodel import PriorDustModel as PhysDustModel

__all__ = ["MB_Model"]


class MB_Model:
    """
    MegaBEAST model that provides member functions to compute
    the likelihood and priors for a specific physical model
    """

    def __init__(self, params):
        self.fd_model = params.fd_model

        # setup the physics model for the beast parameters
        #   uses the same format as the beast priors = megabeast physics model
        # --> needs to be generalized to also handle stellar parameters
        #     define a dict that translates between mb params and physical models
        self.beast_params = ["Av", "Rv", "fA"]
        self.physics_model = {}
        for cparam in self.beast_params:
            cmod = self.fd_model[cparam]
            self.physics_model[cparam] = {"name": cmod["name"]}
            for cname, cval in zip(cmod["varnames"], cmod["varinit"]):
                self.physics_model[cparam][cname] = cval
            # setup the phyics model for this parameter
            self.physics_model[cparam]["model"] = PhysDustModel(
                self.physics_model[cparam]
            )

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
        return (names, values)

    def lnlike(self, phi, star_lnpgriddata, beast_moddata):
        """
        Compute the log(likelihood) for the ensemble parameters

        Parameters
        ----------
        phi: floats
           ensemble parameters

        star_lnpgriddata: dictonary
           contains arrays of the likelihood*grid_weight values and
           indexs to the BEAST model grid

        beast_moddata: dictonary
           contains arrays of the beast parameters for the full beast physics grid

        Returns
        -------
        log(likelihood): float
        """
        # update the values in the physics model using the input ensemble parameters
        k = 0
        cur_physmod = np.full((len(beast_moddata["Av"])), 1.0, dtype=float)
        for cparam in self.beast_params:
            cmod = self.fd_model[cparam]
            if cmod["prior"]["name"] != "fixed":
                for cvar in cmod["varnames"]:
                    self.physics_model[cparam][cvar] = phi[k]
                    k += 1

                # compute the physics model for the full BEAST physics grid
                cur_physmod *= self.physics_model[cparam]["model"](
                    beast_moddata[cparam]
                )

        # compute the each star's integrated probability that it fits the new model
        # including the completeness function
        n_lnps, n_stars = star_lnpgriddata["indxs"].shape
        logintprob = 0.0
        for i in range(n_stars):
            # mask for the finite star's lnpgrid values
            gmask = np.isfinite(star_lnpgriddata["vals"][i])
            # indxs for the star's lnpgrid values in the full beast grid
            curindxs = (star_lnpgriddata["indxs"][i])[gmask]

            # compute the integrated probability
            star_intprob = np.sum(
                (star_lnpgriddata["vals"][i])[gmask]
                * cur_physmod[curindxs]
                * beast_moddata["completeness"][curindxs]
            )
            # checks for spoilers
            if not np.isfinite(star_intprob):
                raise ValueError("invidual integrated star prob is not finite")
            if star_intprob == 0.0:
                raise ValueError("invidual integrated star prob is zero")

            logintprob += np.log(star_intprob)

        return logintprob

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


def lnprob(phi, megabeast_model, star_lnpgriddata, beast_moddata):
    """
    Compute the log(probability) for the ensemble parameters

    Parameters
    ----------
    phi: floats
       ensemble parameters

    megabeast_model : class
        MegaBEAST physical model including priors

    star_lnpgriddata: dictonary
       contains arrays of the likelihood*grid_weight values and
       indexs to the BEAST model grid

    beast_moddata: dictonary
       contains arrays of the beast parameters for the full beast physics grid

    Returns
    -------
    log(probability): float
    """
    ln_prior = megabeast_model.lnprior(phi)

    if not np.isfinite(ln_prior):
        return -np.inf
    return ln_prior + megabeast_model.lnlike(phi, star_lnpgriddata, beast_moddata)


def fit_ensemble(megabeast_model, star_lnpgriddata, beast_moddata):
    """
    Run the MegaBEAST on a single set of BEAST results.

    Parameters
    ----------
    megabeast_model : class
        MegaBEAST physical model including priors

    star_lnpgriddata: dictonary
       contains arrays of the likelihood*grid_weight values and
       indexs to the BEAST model grid

    beast_moddata: dictonary
       contains arrays of the beast parameters for the full beast physics grid

    Returns
    -------
    fit_results : array
        set of best fit parameters
    """
    # standard minimization to find initial values
    def chi2(*args):
        return -1.0 * lnprob(*args)

    result = op.minimize(
        chi2,
        megabeast_model.start_params()[1],
        args=(megabeast_model, star_lnpgriddata, beast_moddata),
        # method="Nelder-Mead",
    )

    # next step would be to
    # run through MCMC to fully sample likelihood
    # maybe include option not to run MCMC

    # print("output")
    # print(megabeast_model.physics_model)
    print(result)

    return result["x"]
