import numpy as np
import argparse
import asdf

from numpy.random import default_rng

from astropy.table import Table, Column

from beast.observationmodel.vega import Vega
from beast.physicsmodel.grid import SEDGrid
from beast.physicsmodel.priormodel import PriorAgeModel, PriorMassModel
from beast.physicsmodel.grid_weights_stars import compute_bin_boundaries
from beast.physicsmodel.grid_and_prior_weights import (
    compute_distance_age_mass_metallicity_weights,
    compute_av_rv_fA_prior_weights,
)
import beast.observationmodel.noisemodel.generic_noisemodel as noisemodel

from megabeast.helpers import read_mbmodel
from megabeast.mbmodel import MBModel

from astropy.table import vstack

__all__ = ["gen_SimObs_from_sedgrid", "simulate_obs"]


def gen_SimObs_from_sedgrid(
    sedgrid,
    sedgrid_noisemodel,
    nsim=0,
    compl_filter="max",
    complcut=None,
    magcut=None,
    ranseed=None,
    vega_fname=None,
    weight_to_use="weight",
    age_prior_model=None,
    mass_prior_model=None,
):
    """
    Generate simulated observations using the physics and observation grids.
    The priors are sampled as they give the ensemble model for the stellar
    and dust distributions (IMF, Av distribution etc.).
    The physics model gives the SEDs based on the priors.
    The observation model gives the noise, bias, and completeness all of
    which are used in simulating the observations.

    Currently written to only work for the toothpick noisemodel.

    Parameters
    ----------
    sedgrid: grid.SEDgrid instance
        model grid

    sedgrid_noisemodel: beast noisemodel instance
        noise model data

    nsim : int
        number of observations to simulate

    compl_filter : str
        Filter to use for completeness (required for toothpick model).
        Set to max to use the max value in all filters.

    complcut : float (defualt=None)
        completeness cut for only including model seds above the cut
        where the completeness cut ranges between 0 and 1.

    magcut : float (defualt=None)
        faint-end magnitude cut for only including model seds brighter
        than the given magnitude in compl_filter.

    ranseed : int
        used to set the seed to make the results reproducable,
        useful for testing

    vega_fname : string
        filename for the vega info, useful for testing

    weight_to_use : string (default='weight')
        Set to either 'weight' (prior+grid), 'prior_weight', 'grid_weight',
        or 'uniform' (this option is valid only when nsim is supplied) to
        choose the weighting for SED selection.

    age_prior_model : dict
        age prior model in the BEAST dictonary format

    mass_prior_model : dict
        mass prior model in the BEAST dictonary format

    Returns
    -------
    simtable : astropy Table
        table giving the simulated observed fluxes as well as the
        physics model parmaeters
    """
    n_models, n_filters = sedgrid.seds.shape
    flux = sedgrid.seds

    # get the vega fluxes for the filters
    _, vega_flux, _ = Vega(source=vega_fname).getFlux(sedgrid.filters)

    # cache the noisemodel values
    model_bias = sedgrid_noisemodel["bias"]
    model_unc = np.fabs(sedgrid_noisemodel["error"])
    model_compl = sedgrid_noisemodel["completeness"]

    # only use models that have non-zero completeness in all filters
    # zero completeness means the observation model is not defined for that filters/flux
    # if complcut is provided, only use models above that completeness cut
    if complcut is not None:
        finalcomplcut = complcut
    else:
        finalcomplcut = 0.0

    ast_defined = model_compl > finalcomplcut
    sum_ast_defined = np.sum(ast_defined, axis=1)
    goodobsmod = sum_ast_defined >= n_filters

    # completeness from toothpick model so n band completeness values
    # require only 1 completeness value for each model
    # max picked to best "simulate" how the photometry detection is done
    if compl_filter.lower() == "max":
        model_compl = np.max(model_compl, axis=1)
    else:
        short_filters = [
            filter.split(sep="_")[-1].upper() for filter in sedgrid.filters
        ]
        if compl_filter.upper() not in short_filters:
            raise NotImplementedError(
                "Requested completeness filter not present:"
                + compl_filter.upper()
                + "\nPossible filters:"
                + "\n".join(short_filters)
            )

        filter_k = short_filters.index(compl_filter.upper())
        print("Completeness from %s" % sedgrid.filters[filter_k])
        model_compl = model_compl[:, filter_k]

    # if magcut is provided, only use models brighter than the magnitude cut
    # in addition to the non-zero completeness criterion
    if magcut is not None:
        fluxcut_compl_filter = 10 ** (-0.4 * magcut) * vega_flux[filter_k]
        goodobsmod = (goodobsmod) & (flux[:, filter_k] >= fluxcut_compl_filter)

    # initialize the random number generator
    rangen = default_rng(ranseed)

    # if the age and mass prior models are given, use them to determine the
    # total number of stars to simulate
    model_indx = np.arange(n_models)
    if (age_prior_model is not None) and (mass_prior_model is not None):
        print("determining the number of simulated stars from age and mass priors")
        nsim = 0
        # logage_range = [min(sedgrid["logA"]), max(sedgrid["logA"])]
        mass_range = [min(sedgrid["M_ini"]), max(sedgrid["M_ini"])]

        # compute the total mass and average mass of a star given the mass_prior_model
        nmass = 100
        masspts = np.logspace(np.log10(mass_range[0]), np.log10(mass_range[1]), nmass)
        mass_prior = PriorMassModel(mass_prior_model)
        massprior = mass_prior(masspts)
        totmass = np.trapz(massprior, masspts)
        avemass = np.trapz(masspts * massprior, masspts) / totmass

        # compute the mass of the remaining stars at each age and
        # simulate the stars assuming everything is complete
        gridweights = sedgrid[weight_to_use]
        gridweights = gridweights / np.sum(gridweights)

        grid_ages = np.unique(sedgrid["logA"])
        age_prior = PriorAgeModel(age_prior_model)
        ageprior = age_prior(grid_ages)
        bin_boundaries = compute_bin_boundaries(grid_ages)
        bin_widths = np.diff(10 ** (bin_boundaries))
        totsim_indx = np.array([], dtype=int)
        for cage, cwidth, cprior in zip(grid_ages, bin_widths, ageprior):
            gmods = sedgrid["logA"] == cage
            cur_mass_range = [
                min(sedgrid["M_ini"][gmods]),
                max(sedgrid["M_ini"][gmods]),
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

        totsimmass = np.sum(sedgrid["M_ini"][totsim_indx])
        print(f"number total simulated stars = {nsim}; mass = {totsimmass}")
        compl_choice = rangen.random(nsim)
        compl_indx = model_compl[totsim_indx] >= compl_choice
        sim_indx = totsim_indx[compl_indx]
        totcompsimmass = np.sum(sedgrid["M_ini"][sim_indx])
        print(
            f"number of simulated stars w/ completeness = {len(sim_indx)}; mass = {totcompsimmass}"
        )

    else:  # total number of stars to simulate set by command line input

        if weight_to_use == "uniform":
            # sample to get the indices of the picked models
            sim_indx = rangen.choice(model_indx[goodobsmod], nsim)

        else:
            gridweights = sedgrid[weight_to_use][goodobsmod] * model_compl[goodobsmod]
            gridweights = gridweights / np.sum(gridweights)

            # sample to get the indexes of the picked models
            sim_indx = rangen.choice(model_indx[goodobsmod], size=nsim, p=gridweights)

        print(f"number of simulated stars = {nsim}")

    # setup the output table
    ot = Table()
    qnames = list(sedgrid.keys())
    # simulated data
    for k, filter in enumerate(sedgrid.filters):
        simflux_wbias = flux[sim_indx, k] + model_bias[sim_indx, k]

        simflux = rangen.normal(loc=simflux_wbias, scale=model_unc[sim_indx, k])

        bname = filter.split(sep="_")[-1].upper()
        fluxname = f"{bname}_FLUX"
        colname = f"{bname}_RATE"
        magname = f"{bname}_VEGA"
        ot[fluxname] = Column(simflux)
        ot[colname] = Column(ot[fluxname] / vega_flux[k])
        pindxs = ot[colname] > 0.0
        nindxs = ot[colname] <= 0.0
        ot[magname] = Column(ot[colname])
        ot[magname][pindxs] = -2.5 * np.log10(ot[colname][pindxs])
        ot[magname][nindxs] = 99.999

        # add in the physical model values in a form similar to
        # the output simulated (physics+obs models) values
        # useful if using the simulated data to interpolate ASTs
        #   (e.g. for MATCH)
        fluxname = f"{bname}_INPUT_FLUX"
        ratename = f"{bname}_INPUT_RATE"
        magname = f"{bname}_INPUT_VEGA"
        ot[fluxname] = Column(flux[sim_indx, k])
        ot[ratename] = Column(ot[fluxname] / vega_flux[k])
        pindxs = ot[ratename] > 0.0
        nindxs = ot[ratename] <= 0.0
        ot[magname] = Column(ot[ratename])
        ot[magname][pindxs] = -2.5 * np.log10(ot[ratename][pindxs])
        ot[magname][nindxs] = 99.999

    # model parmaeters
    for qname in qnames:
        ot[qname] = Column(sedgrid[qname][sim_indx])

    return ot


def simulate_obs(
    physgrid_list,
    noise_model_list,
    output_catalog=None,
    beastinfo_list=None,
    mbmodel=None,
    nsim=0,
    compl_filter="max",
    complcut=None,
    magcut=None,
    weight_to_use="weight",
    ranseed=None,
):
    """
    Simulate photometry based on a physicsmodel grid(s) and observation model(s).

    Parameters
    ----------
    physgrid_list : list of strings
        Name of the physics model file.  If there are multiple physics model
        grids (i.e., if there are subgrids), list them all here, and they will
        each be sampled nsim/len(physgrid_list) times.

    noise_model_list : list of strings
        Name of the noise model file.  If there are multiple files for
        physgrid_list (because of subgrids), list the noise model files
        associated with each physics model file.

    output_catalog : string
        Name of the output simulated photometry catalog

    beastinfo_list : list of strings (default=None)
        Name of the beast info file.  The mass and age prior models are read
        from this model to use to compute the number of stars to simulate. If
        there are multiple files for physgrid_list (because of subgrids), list
        the beast info files associated with each physics model file.
        Cannot be used at the same time as mbmodel.

    mbmodel: megabeast.MBModel
        MegaBEAST model that provides a full ensemble physics model.
        Cannot be used at the same time as the beastinfo_list.

    n_sim : int (default=100)
        Number of simulated objects to create if beastinfo_list is not given. If
        nsim/len(physgrid_list) isn't an integer, this will be increased so that
        each grid has the same number of samples.

    compl_filter : str (default=max)
        filter to use for completeness (required for toothpick model)
        set to max to use the max value in all filters

    complcut : float (defualt=None)
        completeness cut for only including model seds above the given
        completeness, where the cut ranges from 0 to 1.

    magcut : float (defualt=None)
        faint-end magnitude cut for only including model seds brighter
        than the given magnitude in compl_filter.

    weight_to_use : string (default='weight')
        Set to either 'weight' (prior+grid), 'prior_weight', 'grid_weight',
        or 'uniform' (this option is valid only when nsim is supplied) to
        choose the weighting for SED selection.

    ranseed : int
        seed for random number generator
    """
    if (beastinfo_list is not None) & (mbmodel is not None):
        raise Exception("beastinfo_list and mbmodel cannot both be set")

    # numbers of samples to do
    # (ensure there are enough for even sampling of multiple model grids)
    n_phys = len(np.atleast_1d(physgrid_list))
    nsim = int(nsim)
    samples_per_grid = int(np.ceil(nsim / n_phys))

    if complcut is not None:
        complcut = float(complcut)

    if magcut is not None:
        magcut = float(magcut)

    if ranseed is not None:
        ranseed = int(ranseed)

    # list to hold all simulation tables
    simtable_list = []

    # make a table for each physics model + noise model
    for k, physgrid in enumerate(np.atleast_1d(physgrid_list)):
        noise_model = np.atleast_1d(noise_model_list)[k]

        # get the physics model grid - includes priors
        modelsedgrid = SEDGrid(str(physgrid))

        # read in the noise model - includes bias, unc, and completeness
        noisegrid = noisemodel.get_noisemodelcat(str(noise_model))

        if beastinfo_list is not None:
            with asdf.open(np.atleast_1d(beastinfo_list)[k]) as af:
                binfo = af.tree
                age_prior_model = binfo["age_prior_model"]
                mass_prior_model = binfo["mass_prior_model"]
        else:
            age_prior_model = None
            mass_prior_model = None

        # update the prior_weights if mbmodel is set
        if mbmodel is not None:

            if nsim > 0:
                age_prior_model = None
                mass_prior_model = None
            else:
                age_prior_model = (mbmodel.physics_model["logA"],)
                mass_prior_model = (mbmodel.physics_model["M_ini"],)

            # mock up a beast spectral grid so that the BEAST code can be used
            cur_physmod = Table()
            for cparam in ["logA", "M_ini", "Z", "distance", "Av", "Rv", "f_A"]:
                cur_physmod[cparam] = modelsedgrid.grid[cparam]
            n_mods = len(cur_physmod["logA"])
            for cparam in ["grid_weight", "prior_weight", "weight"]:
                cur_physmod[cparam] = np.full(n_mods, 1.0)

            print("computing megabeast ensemble physics model (replaces beast priors)")

            dust_prior = compute_av_rv_fA_prior_weights(
                cur_physmod["Av"],
                cur_physmod["Rv"],
                cur_physmod["f_A"],
                cur_physmod["distance"],
                av_prior_model=mbmodel.physics_model["Av"],
                rv_prior_model=mbmodel.physics_model["Rv"],
                fA_prior_model=mbmodel.physics_model["f_A"],
            )

            compute_distance_age_mass_metallicity_weights(
                cur_physmod,
                distance_prior_model={"name": "flat"},
                age_prior_model=mbmodel.physics_model["logA"]["model"],
                mass_prior_model=mbmodel.physics_model["M_ini"]["model"],
                met_prior_model={"name": "flat"},
            )

            cur_physmod["prior_weight"] *= dust_prior

            modelsedgrid.grid["prior_weight"] = cur_physmod["prior_weight"].data
            modelsedgrid.grid["grid_weight"] = cur_physmod["grid_weight"].data
            modelsedgrid.grid["weight"] = (
                modelsedgrid.grid["prior_weight"] * modelsedgrid.grid["grid_weight"]
            )

        # generate the table
        simtable = gen_SimObs_from_sedgrid(
            modelsedgrid,
            noisegrid,
            age_prior_model=age_prior_model,
            mass_prior_model=mass_prior_model,
            nsim=samples_per_grid,
            compl_filter=compl_filter,
            complcut=complcut,
            magcut=magcut,
            weight_to_use=weight_to_use,
            ranseed=ranseed,
        )

        # append to the list
        simtable_list.append(simtable)

    # stack all the tables into one and write it out
    fincat = vstack(simtable_list)
    if output_catalog is not None:
        fincat.write(output_catalog, overwrite=True)

    return fincat


def main():
    # commandline parser
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--physgrid_list",
        "-p",
        metavar="PHYS_MODEL",
        required=True,
        nargs="+",
        help="filename(s) of physics grid(s)",
    )
    parser.add_argument(
        "--noise_model_list",
        "-n",
        metavar="NOISE_MODEL",
        required=True,
        nargs="+",
        help="filename(s) of observation/noise grid(s)",
    )
    parser.add_argument(
        "--output_catalog",
        "-c",
        required=True,
        help="filename for simulated observations",
    )
    parser.add_argument(
        "--beastinfo_list",
        metavar="BEAST_INFO",
        nargs="+",
        help="filename(s) of beast info file(s)",
    )
    parser.add_argument("--ensembleparams", help="ensemble parameters file")
    parser.add_argument(
        "--nsim", default=0, type=int, help="number of simulated objects"
    )
    parser.add_argument(
        "--compl_filter",
        default="F475W",
        help="filter for completeness, set to max for max of values from all filters",
    )
    parser.add_argument(
        "--complcut",
        default=None,
        type=float,
        help="completeness cut for selecting seds above the completeness cut",
    )
    parser.add_argument(
        "--magcut",
        default=None,
        type=float,
        help="magnitdue cut for selecting seds brighter than the magcut in compl_filter",
    )
    parser.add_argument(
        "--weight_to_use",
        default="weight",
        type=str,
        help="""Set to either 'weight' (prior+grid), 'prior_weight', or
        'grid_weight' to choose the weighting for SED selection.""",
    )
    parser.add_argument(
        "--ranseed", default=None, type=int, help="seed for random number generator"
    )
    args = parser.parse_args()

    # read in the megabeast model if supplied
    if args.mbmodel:
        params = read_mbmodel(args.mbmodel)
        mbmodel = MBModel(params.stellar_model, params.dust_model)
    else:
        mbmodel = None

    # run observation simulator
    simulate_obs(
        args.physgrid_list,
        args.noise_model_list,
        output_catalog=args.output_catalog,
        beastinfo_list=args.beastinfo_list,
        mbmodel=mbmodel,
        nsim=args.nsim,
        compl_filter=args.compl_filter,
        complcut=args.complcut,
        magcut=args.magcut,
        weight_to_use=args.weight_to_use,
        ranseed=args.ranseed,
    )


if __name__ == "__main__":  # pragma: no cover

    main()
