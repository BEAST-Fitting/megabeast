import argparse
import h5py
import numpy as np
import asdf

from beast.physicsmodel.grid import SEDGrid

from megabeast.mbsettings import mbsettings
from megabeast.singlepop_dust_model import MB_Model, fit_ensemble


def main():
    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "settings_file",
        help="Name of the file that contains the MegaBEAST settings",
    )
    args = parser.parse_args()

    # read in the parameters
    mbparams = mbsettings(args.settings_file)
    print(mbparams.fd_model)
    exit()

    sedsfile = mbparams.beast_base + "_seds.grid.hd5"
    obsmodfile = mbparams.beast_base + "_noisemodel.grid.hd5"
    infofile = mbparams.beast_base + "_beast_info.asdf"

    lnpfile = mbparams.beast_base + "_lnp.hd5"

    # get the beast data needed
    #  using SEDGrid as it is faster than beast.tools.read_beast_data.read_sed_data
    beast_data = {}
    sgrid = SEDGrid(sedsfile, backend="disk")
    for cparam in mbparams.param_list:
        beast_data[cparam] = sgrid.grid[cparam]

    # use the maximum completeness across the bands as the correct obsmodel
    # would only have one completeness value per model
    # max is the best approximation for the toothpick model
    with h5py.File(obsmodfile, "r") as obs_hdf:
        beast_data["completeness"] = np.max(obs_hdf["completeness"], axis=1)

    # BEAST prior model from the beast_info file
    beast_pmodel = {}
    with asdf.open(infofile) as af:
        # tree = af.tree
        beast_pmodel["AV"] = af.tree["av_prior_model"]
        beast_pmodel["RV"] = af.tree["rv_prior_model"]
        beast_pmodel["fA"] = af.tree["fA_prior_model"]

    print(beast_pmodel)
    exit()

    megabeast_model = MB_Model(mbparams)

    bestparams = fit_ensemble(beast_data, lnpfile, beast_pmodel, megabeast_model)
    print(bestparams)


if __name__ == "__main__":

    main()
