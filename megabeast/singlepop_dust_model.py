import numpy as np

# import scipy.optimize as op
import emcee
from numpy.random import default_rng
import scipy

from beast.physicsmodel.priormodel import PriorAgeModel as PhysAgeModel
from beast.physicsmodel.priormodel import PriorMassModel as PhysMassModel
from beast.physicsmodel.priormodel import PriorDustModel as PhysDustModel
from beast.physicsmodel.grid_weights_stars import compute_bin_boundaries

__all__ = ["MB_Model"]


def _get_predicted_num_stars(bphysparams, bphysmod, physmodage, physmodmass):
    """
    Calculate the expected number of stars based on the physics model as
    given on the BEAST model grid including completeness.

    Parameters
    ----------
    bphysparams : astropy.table
        table giving the physical parameters, weights, and completeness on
        the BEAST physical grid

    bphysmod : array
        probability of the physics model on the BEAST physics grid

    physmodage : beast.physicmodel.priormodel
        physics model of for age

    physmodmass : beast.physicmodel.priormodel
        physics model of for initial mass

    Returns
    -------
    n_stars : float
        number of stars expected
    """

    # initialize the random number generator
    rangen = default_rng()

    model_indx = np.arange(len(bphysparams["M_ini"]))

    nsim = 0
    # logage_range = [min(sedgrid["logA"]), max(sedgrid["logA"])]
    mass_range = [min(bphysparams["M_ini"]), max(bphysparams["M_ini"])]

    # compute the total mass and average mass of a star given the mass_prior_model
    nmass = 100
    masspts = np.logspace(np.log10(mass_range[0]), np.log10(mass_range[1]), nmass)
    massprior = physmodmass(masspts)
    totmass = np.trapz(massprior, masspts)
    avemass = np.trapz(masspts * massprior, masspts) / totmass

    # compute the mass of the remaining stars at each age and
    # simulate the stars assuming everything is complete
    gridweights = bphysmod * bphysparams["grid_weight"]
    gridweights = gridweights / np.sum(gridweights)

    grid_ages = np.unique(bphysparams["logA"])
    ageprior = physmodage(grid_ages)
    bin_boundaries = compute_bin_boundaries(grid_ages)
    bin_widths = np.diff(10 ** (bin_boundaries))
    totsim_indx = np.array([], dtype=int)
    for cage, cwidth, cprior in zip(grid_ages, bin_widths, ageprior):
        gmods = bphysparams["logA"] == cage
        cur_mass_range = [
            min(bphysparams["M_ini"][gmods]),
            max(bphysparams["M_ini"][gmods]),
        ]
        gmass = (masspts >= cur_mass_range[0]) & (masspts <= cur_mass_range[1])
        curmasspts = masspts[gmass]
        curmassprior = massprior[gmass]
        totcurmass = np.trapz(curmassprior, curmasspts)

        # compute the mass remaining at each age -> this is the mass to simulate
        simmass = cprior * cwidth * totcurmass / totmass
        nsim_curage = int(round(simmass / avemass))

        # simluate the stars at the current age
        curweights = gridweights[gmods]
        if np.sum(curweights) > 0:
            curweights /= np.sum(curweights)
            cursim_indx = rangen.choice(
                model_indx[gmods], size=nsim_curage, p=curweights
            )

            totsim_indx = np.concatenate((totsim_indx, cursim_indx))

            nsim += nsim_curage
        # totsimcurmass = np.sum(sedgrid["M_ini"][cursim_indx])
        # print(cage, totcurmass / totmass, simmass, totsimcurmass, nsim_curage)

    # totsimmass = np.sum(bphysparams["M_ini"][totsim_indx])
    # print(f"number total simulated stars = {nsim}; mass = {totsimmass}")
    compl_choice = rangen.random(nsim)
    compl_indx = bphysparams["completeness"][totsim_indx] >= compl_choice
    sim_indx = totsim_indx[compl_indx]
    # totcompsimmass = np.sum(bphysparams["M_ini"][sim_indx])
    # print(
    #     f"number of simulated stars w/ completeness = {len(sim_indx)}; mass = {totcompsimmass}"
    # )

    return len(sim_indx)


class MB_Model:
    """
    MegaBEAST model that provides member functions to compute
    the likelihood and priors for a specific physical model
    """

    def __init__(self, params):
        self.star_model = params.stellar_model
        self.dust_model = params.fd_model

        # setup the physics model for the beast parameters
        #   uses the same format as the beast priors = megabeast physics model
        # --> needs to be generalized to also handle stellar parameters
        #     define a dict that translates between mb params and physical models
        self.params = ["logA", "M_ini", "Av", "Rv", "fA"]
        self.physics_model = {}
        for cparam in self.params:
            if cparam in self.star_model.keys():
                cmod = self.star_model[cparam]
            elif cparam in self.dust_model.keys():
                cmod = self.dust_model[cparam]
            else:
                raise ValueError("requested parameter not in mbsetting file")

            self.physics_model[cparam] = {"name": cmod["name"]}
            self.physics_model[cparam]["varnames"] = cmod["varnames"]
            self.physics_model[cparam]["prior"] = cmod["prior"]
            for cname, cval in zip(cmod["varnames"], cmod["varinit"]):
                self.physics_model[cparam][cname] = cval

            # setup the phyics model for this parameter
            if cparam in self.star_model.keys():
                if cparam == "logA":
                    self.physics_model[cparam]["x"] = self.star_model["x"]
                    self.physics_model[cparam]["nsubvars"] = len(
                        self.physics_model[cparam]["x"]
                    )
                    self.physics_model[cparam]["model"] = PhysAgeModel(
                        self.physics_model[cparam]
                    )
                elif cparam == "M_ini":
                    # currently no parameters possible
                    self.physics_model[cparam]["model"] = PhysMassModel(
                        self.physics_model[cparam]
                    )
            elif cparam in self.dust_model.keys():
                self.physics_model[cparam]["nsubvars"] = 1
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
        mod = self.physics_model
        names = []
        values = []
        for ckey in mod.keys():
            if mod[ckey]["prior"]["name"] != "fixed":
                for k, cparam in enumerate(mod[ckey]["varnames"]):
                    if (
                        len(np.atleast_1d(mod[ckey][cparam])) > 1
                    ):  # expand into multiple parameters
                        for l, cval in enumerate(mod[ckey][cparam]):
                            names.append(f"{ckey}_{cparam}{l+1}")
                            values.append(cval)
                    else:
                        names.append(f"{ckey}_{cparam}")
                        values.append(mod[ckey][cparam])

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
        for cparam in self.params:
            cmod = self.physics_model[cparam]
            if cmod["prior"]["name"] != "fixed":
                for cvar in cmod["varnames"]:
                    if cmod["nsubvars"] > 1:
                        for j in range(cmod["nsubvars"]):
                            self.physics_model[cparam]["values"][j] = phi[k]
                            k += 1
                    else:
                        self.physics_model[cparam][cvar] = phi[k]
                        k += 1

                # compute the physics model for the full BEAST physics grid
                cur_physmod *= self.physics_model[cparam]["model"](
                    beast_moddata[cparam]
                )
                # if cparam == "logA":
                #     print(self.physics_model[cparam]["values"])

        n_lnps, n_stars = star_lnpgriddata["indxs"].shape

        # compute probability the model produces the observed number of stars
        pred_stars = _get_predicted_num_stars(
            beast_moddata,
            cur_physmod,
            self.physics_model["logA"]["model"],
            self.physics_model["M_ini"]["model"],
        )
        # cacluate the probability of the observed number of stars
        #  ln form based on equation 8 in Weisz et al. (2013, ApJ, 762, 123)
        logintprob = (
            n_stars * np.log(pred_stars)
            - pred_stars
            - scipy.special.gammaln(n_stars + 1)
        )
        # print(pred_stars, n_stars, logintprob)

        # compute the each star's integrated probability that it fits the new model
        # including the completeness function
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
        cmod = self.physics_model
        k = 0
        for cparam in cmod.keys():
            cprior = cmod[cparam]["prior"]
            cname = cprior["name"]
            if cname != "fixed":
                if cname == "flat":
                    for i, vparam in enumerate(cmod[cparam]["varnames"]):
                        if cmod[cparam]["nsubvars"] > 1:
                            for j in range(len(np.atleast_1d(cprior["minmax"][i]))):
                                if (
                                    phi[k] < cprior["minmax"][0][j]
                                    or phi[k] > cprior["minmax"][1][j]
                                ):
                                    return -np.inf
                                k += 1
                        else:
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


def _get_best_fit_params(sampler):
    """
    Determine the best fit parameters given an emcee sampler object
    """
    # very likely a faster way
    max_lnp = -1e6
    nwalkers, nsteps = sampler.lnprobability.shape
    for k in range(nwalkers):
        tmax_lnp = np.nanmax(sampler.lnprobability[k])
        if tmax_lnp > max_lnp:
            max_lnp = tmax_lnp
            (indxs,) = np.where(sampler.lnprobability[k] == tmax_lnp)
            fit_params_best = sampler.chain[k, indxs[0], :]

    return fit_params_best


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

    sparams = megabeast_model.start_params()[1]

    # result = op.minimize(
    # result = op.least_squares(
    #    chi2,
    #    sparams,
    #    args=(megabeast_model, star_lnpgriddata, beast_moddata),
    #    ftol=1e-20,
    #    xtol=1e-20
    #    method="Nelder-Mead",
    # )

    ndim, nwalkers = len(sparams), 5 * len(sparams)

    pos = sparams * (1.0 + 1e-1 * np.random.randn(nwalkers, ndim))

    sampler = emcee.EnsembleSampler(
        nwalkers, ndim, lnprob, args=(megabeast_model, star_lnpgriddata, beast_moddata)
    )
    nsteps = 100
    sampler.run_mcmc(pos, nsteps, progress=True)
    # samples = sampler.get_chain()

    # next step would be to
    # run through MCMC to fully sample likelihood
    # maybe include option not to run MCMC

    # print("output")
    # print(megabeast_model.physics_model)
    # print(result)

    return _get_best_fit_params(sampler)
    # return result["x"]
