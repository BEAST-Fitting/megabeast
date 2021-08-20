import numpy as np
from numpy.random import default_rng

from beast.tools.read_beast_data import read_lnp_data
from beast.physicsmodel.grid_weights_stars import compute_bin_boundaries

__all__ = [
    "get_likelihoods",
    "precompute_mass_multipliers",
    "get_predicted_num_stars",
    "get_predicted_num_stars_simulate",
]


def get_likelihoods(ppdf_file, beast_model_data):
    """
    Read in the saved BEAST sparse posterior PDFs and divide by the BEAST
    priors to get the spare liklilhoods.

    Parameters
    ----------
    ppdf_file : string
        filename of the saved BEAST posterior PDFs

    Returns
    -------
    lnp_data : dictonary
       contains arrays of the likelihood values (vals) and indices to the BEAST model grid (indxs)
    """
    # BEAST saves posterior PDFs labeled as log(pPDF)
    lnpdata = read_lnp_data(ppdf_file)

    # divide by the BEAST prior weights to return to recover the likelihoods
    n_lnps, n_stars = lnpdata["indxs"].shape
    for i in range(n_stars):
        indxs = lnpdata["indxs"][:, i]
        gmask = np.isfinite(lnpdata["vals"][:, i])
        lnpdata["vals"][gmask, i] = (
            np.exp(lnpdata["vals"][gmask, i])
            / beast_model_data["prior_weight"][indxs[gmask]]
        )

    return lnpdata


def precompute_mass_multipliers(bphysparams, physmodmass):
    """
    Calculate the value to mulitply the SFR to get the total mass in stars
    at all ages and masses given the physics model on the BEAST model grid.

    Parameters
    ----------
    bphysparams : astropy.table
        table giving the physical parameters, weights, and completeness on
        the BEAST physical grid

    physmodmass : beast.physicmodel.priormodel
        physics model of for initial mass

    Returns
    -------
    sfhmassinfo : dict
        "massmult" gives the value to muliply the SFH at each age and metallicity
        "ages" gives the ages and "Zs" gives the metallicities
    """
    mass_range = [min(bphysparams["M_ini"]), max(bphysparams["M_ini"])]

    # compute the total mass and average mass of a star given the mass_prior_model
    nmass = 100
    masspts = np.logspace(np.log10(mass_range[0]), np.log10(mass_range[1]), nmass)
    massprior = physmodmass(masspts)
    totmass = np.trapz(massprior, masspts)
    avemass = np.trapz(masspts * massprior, masspts) / totmass

    # loop over all the ages and compute the mass to simulate
    #   ***need to add metallicity as well***
    grid_ages = np.unique(bphysparams["logA"])
    bin_boundaries = compute_bin_boundaries(grid_ages)
    bin_widths = np.diff(10 ** (bin_boundaries))
    grid_mets = np.unique(bphysparams["Z"])
    massmults = np.full((len(grid_ages), len(grid_mets)), 0.0)
    for i, cage in enumerate(grid_ages):
        for j, cmet in enumerate(grid_mets):
            gmods = (bphysparams["logA"] == cage) & (bphysparams["Z"] == cmet)
            cur_mass_range = [
                min(bphysparams["M_ini"][gmods]),
                max(bphysparams["M_ini"][gmods]),
            ]
            gmass = (masspts >= cur_mass_range[0]) & (masspts <= cur_mass_range[1])
            curmasspts = masspts[gmass]
            curmassprior = massprior[gmass]
            totcurmass = np.trapz(curmassprior, curmasspts)

            # compute the mass remaining at each age -> this is the mass to simulate
            massmults[i, j] = bin_widths[i] * totcurmass / totmass

    return {
        "massmult": massmults,
        "ages": grid_ages,
        "Zs": grid_mets,
        "avemass": avemass,
    }


def get_predicted_num_stars(
    massmult_info, bphysparams, bphysmod, physmodage
):
    """
    Calculate the expected number of stars based on the physics model as
    given on the BEAST model grid including completeness.

    Parameters
    ----------
    massmult_info : dict
        "massmult" gives the value to muliply the SFH at each age and metallicity
        "ages" gives the ages and "Zs" gives the metallicities

    bphysparams : astropy.table
        table giving the physical parameters, weights, and completeness on
        the BEAST physical grid

    bphysmod : array
        probability of the physics model on the BEAST physics grid

    physmodage : beast.physicmodel.priormodel
        physics model of for age

    Returns
    -------
    n_stars : float
        number of stars expected
    """
    gridweights = bphysmod * bphysparams["grid_weight"]
    gridweights = gridweights / np.sum(gridweights)
    ageprior = physmodage(massmult_info["ages"])
    # assumes a flat metallicity prior -> will revise when metallicity
    metprior = np.full(len(massmult_info["Zs"]), 1.0)
    metprior /= np.sum(metprior)

    # loop over the age and mass points and compute the number of stars expected
    # at each point including the completeness
    n_totstars = 0.0
    for i, cage in enumerate(massmult_info["ages"]):
        for j, cmet in enumerate(massmult_info["Zs"]):
            gmods = (bphysparams["logA"] == cage) & (bphysparams["Z"] == cmet)

            # compute the mass remaining at each age -> this is the mass to simulate
            simmass = ageprior[i] * metprior[j] * massmult_info["massmult"][i, j]

            # number of stars (can be fractional)
            nstars_curage = simmass / massmult_info["avemass"]

            # normalize the total probability to have number of stars
            curweights = gridweights[gmods]

            if np.sum(curweights) > 0:
                curweights *= nstars_curage / np.sum(curweights)
                # account for the completeness
                curweights *= bphysparams["completeness"][gmods]
                nstars_curage_wcomp = np.sum(curweights)

                n_totstars += nstars_curage_wcomp

    return n_totstars


def get_predicted_num_stars_simulate(
    massmult_info, bphysparams, bphysmod, physmodage
):
    """
    Calculate the expected number of stars based on the physics model as
    given on the BEAST model grid including completeness.

    Parameters
    ----------
    massmult_info : dict
        "massmult" gives the value to muliply the SFH at each age and metallicity
        "ages" gives the ages and "Zs" gives the metallicities

    bphysparams : astropy.table
        table giving the physical parameters, weights, and completeness on
        the BEAST physical grid

    bphysmod : array
        probability of the physics model on the BEAST physics grid

    physmodage : beast.physicmodel.priormodel
        physics model of for age

    Returns
    -------
    n_stars : float
        number of stars expected
    """
    # initialize the random number generator
    rangen = default_rng()

    # setup
    model_indx = np.arange(len(bphysparams["M_ini"]))
    nsim = 0

    gridweights = bphysmod * bphysparams["grid_weight"]
    gridweights = gridweights / np.sum(gridweights)
    ageprior = physmodage(massmult_info["ages"])
    # assumes a flat metallicity prior -> will revise when metallicity
    metprior = np.full(len(massmult_info["Zs"]), 1.0)
    metprior /= np.sum(metprior)

    # loop over the age and mass points and compute the number of stars expected
    # at each point including the completeness
    totsim_indx = np.array([], dtype=int)
    for i, cage in enumerate(massmult_info["ages"]):
        for j, cmet in enumerate(massmult_info["Zs"]):
            gmods = (bphysparams["logA"] == cage) & (bphysparams["Z"] == cmet)

            # compute the mass remaining at each age -> this is the mass to simulate
            simmass = ageprior[i] * metprior[j] * massmult_info["massmult"][i, j]
            nsim_curage = int(simmass / massmult_info["avemass"])

            # simluate the stars at the current age and metallicity
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
