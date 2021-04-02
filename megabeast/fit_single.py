import argparse
import h5py
import numpy as np

from beast.physicsmodel.grid import SEDGrid

from megabeast.mbsettings import mbsettings
from megabeast.helpers import get_likelihoods
from megabeast.singlepop_dust_model import MB_Model, fit_ensemble


def main():
    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "settings_file",
        help="Name of the file that contains the MegaBEAST settings",
    )
    args = parser.parse_args()

    # read the parameters into a class
    mbparams = mbsettings(args.settings_file)
    # should validate the inputs - make sure the required info is provided

    # BEAST files used by MegaBEAST
    sedsfile = mbparams.beast_base + "_seds.grid.hd5"
    obsmodfile = mbparams.beast_base + "_noisemodel.grid.hd5"
    lnpfile = mbparams.beast_base + "_lnp.hd5"

    # get the BEAST physics model info needed
    #  using SEDGrid as it is faster than beast.tools.read_beast_data.read_sed_data
    #  only read in the columns specifically needed
    beast_mod_data = {}
    beast_physmod_param_list = ["Av", "Rv", "f_A", "M_ini", "logA", "Z", "distance", "prior_weight"]
    sgrid = SEDGrid(sedsfile, backend="disk")
    for cparam in beast_physmod_param_list:
        beast_mod_data[cparam] = sgrid.grid[cparam]

    # get the completness from the BEAST observation model
    #   use the maximum completeness across the bands as the correct obsmodel
    #   would only have one completeness value per model
    #   max is the best approximation for the toothpick model (maybe???  average??)
    with h5py.File(obsmodfile, "r") as obs_hdf:
        beast_mod_data["completeness"] = np.max(obs_hdf["completeness"], axis=1)

    # get the saved nD (spare) posterior PDF for each star
    lnp_data = get_likelihoods(lnpfile, beast_mod_data)

    exit()

    megabeast_model = MB_Model(mbparams)

    bestparams = fit_ensemble(beast_data, lnpfile, megabeast_model)
    print(bestparams)


if __name__ == "__main__":

    main()
