# system

# other packages
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import fits

# beast
from beast.physicsmodel.prior_weights_dust import PriorWeightsDust
from beast.tools import read_beast_data

# megabeast
from .read_megabeast_input import read_megabeast_input
from .ensemble_model import lnprob, _two_lognorm

import pdb


def plot_input_data(megabeast_input_file, chi2_plot=[], log_scale=False):
    """
    Parameters
    ----------
    megabeast_input_file : string
        Name of the file that contains settings, filenames, etc

    chi2_plot : list of floats (default=[])
        Make A_V histogram(s) with chi2 less than each of the values in this list

    log_scale : boolean (default=False)
        If True, make the histogram x-axis a log scale (to visualize log-normal
        A_V distribution)

    """

    # read in the settings from the file
    mb_settings = read_megabeast_input(megabeast_input_file)

    # get the project name
    projectname = mb_settings["projectname"]

    # read in the beast data that is needed by all the pixels
    beast_data = {}
    # - SED data
    beast_data.update(read_beast_data.read_sed_data(
        mb_settings["beast_seds_filename"],
        param_list=["Av"]#, "Rv", "f_A"]
    ))
    # - max completeness
    beast_data.update(read_beast_data.read_noise_data(
        mb_settings["beast_noise_filename"],
        param_list=["completeness"],
    ))
    beast_data["completeness"] = np.max(beast_data["completeness"], axis=1)

    # read in the nstars image
    nstars_image, nstars_header = fits.getdata(
        mb_settings["nstars_filename"], header=True
    )
    # dimensions of images/plotting
    y_dimen = nstars_image.shape[0]
    x_dimen = nstars_image.shape[1]




    # set up multi-page figure
    if not log_scale:
        pp = PdfPages("{0}_megabeast/plot_input_data.pdf".format(projectname))
    if log_scale:
        pp = PdfPages("{0}_megabeast/plot_input_data_log.pdf".format(projectname))

    # save the best-fit A_V
    best_av = [[[] for j in range(x_dimen)] for i in range(y_dimen)]
    best_av_chi2 = [[[] for j in range(x_dimen)] for i in range(y_dimen)]

    # -----------------
    # Completeness vs A_V
    # -----------------

    print("")
    print("Making completeness/Av plot")
    print("")

    # set up figure
    fig = plt.figure(figsize=(6, 6))
    plt.subplot(1, 1, 1)

    for i in tqdm(range(y_dimen), desc="y pixels"):
        for j in tqdm(range(x_dimen), desc="x pixels"):
            # for i in tqdm(range(int(y_dimen/3)), desc='y pixels'):
            #    for j in tqdm(range(int(x_dimen/3)), desc='x pixels'):
            # for i in [0]:
            #    for j in [12]:

            if nstars_image[i, j] > 20:

                # get info about the fits
                lnp_filename = mb_settings["lnp_file_prefix"]+"_{0}_{1}_lnp.hd5".format(j, i)
                lnp_data = read_beast_data.read_lnp_data(
                    lnp_filename,
                    nstars=nstars_image[i,j],
                    shift_lnp=True,
                )

                # get the completeness and BEAST model parameters for the
                #   same grid points as the sparse likelihoods
                lnp_grid_vals = read_beast_data.get_lnp_grid_vals(beast_data, lnp_data)

                # grab the things we want to plot
                plot_av = lnp_grid_vals["Av"]
                plot_comp = lnp_grid_vals["completeness"]

                for n in range(nstars_image[i,j]):

                    # plot a random subset of the AVs and completenesses
                    if (i % 3 == 0) and (j % 3 == 0):
                        plot_these = np.random.choice(
                            plot_av[:, n].size, size=20, replace=False
                        )
                        plt.plot(
                            plot_av[plot_these, n]
                            + np.random.normal(scale=0.02, size=plot_these.size),
                            plot_comp[plot_these, n],
                            marker=".",
                            c="black",
                            ms=3,
                            mew=0,
                            linestyle="None",
                            alpha=0.05,
                        )

                    # also overplot the values for the best fit
                    max_ind = np.where(
                        lnp_data["vals"][:, n] == np.max(lnp_data["vals"][:, n])
                    )[0][0]
                    best_av[i][j].append(plot_av[max_ind, n])
                    best_av_chi2[i][j].append(-2 * np.max(lnp_data["vals"][:, n]))
                    if (i % 3 == 0) and (j % 3 == 0):
                        plt.plot(
                            plot_av[max_ind, n] + np.random.normal(scale=0.01),
                            plot_comp[max_ind, n],
                            marker=".",
                            c="magenta",
                            ms=2,
                            mew=0,
                            linestyle="None",
                            alpha=0.3,
                            zorder=9999,
                        )
                    # pdb.set_trace()

    ax = plt.gca()
    ax.set_xlabel(r"$A_V$")
    ax.set_ylabel("Completeness")

    pp.savefig()

    # -----------------
    # histograms of AVs
    # -----------------

    print("")
    print("Making Av Histograms")
    print("")

    # set up figure
    fig = plt.figure(figsize=(x_dimen * 2, y_dimen * 2))

    # flat list of A_V
    # https://stackoverflow.com/questions/952914/making-a-flat-list-out-of-list-of-lists-in-python
    flat_av = [i for sublist in best_av for item in sublist for i in item]
    # grab the max A_V of all of them
    # max_av = max(flat_av)
    # define bins
    if not log_scale:
        uniq_av = np.unique(flat_av)
        gap = np.min(np.diff(uniq_av))
        bins = np.arange(uniq_av[0], uniq_av[-1], gap)
    if log_scale:
        uniq_av = np.unique(np.log10(flat_av))
        gap = (uniq_av[-1] - uniq_av[0]) / len(uniq_av)
        bins = np.arange(uniq_av[0], uniq_av[-1], gap)

    for i in tqdm(range(y_dimen), desc="y pixels"):
        for j in tqdm(range(x_dimen), desc="x pixels"):
    # for i in [0]:
    #    for j in [12]:

            if nstars_image[i, j] > 20:

                # set up the subplot
                plt.subplot(y_dimen, x_dimen, (y_dimen - i - 1) * (x_dimen) + j + 1)

                # make a histogram
                if best_av[i][j] != []:
                    if not log_scale:
                        h = plt.hist(
                            best_av[i][j],
                            bins=bins.size,
                            range=(uniq_av[0] - gap / 2, uniq_av[-1] + gap / 2),
                            facecolor="xkcd:azure",
                            linewidth=0.25,
                            edgecolor="xkcd:azure",
                        )
                    if log_scale:
                        h = plt.hist(
                            np.log10(best_av[i][j]),
                            bins=bins.size,
                            range=(uniq_av[0] - gap / 2, uniq_av[-1] + gap / 2),
                            facecolor="xkcd:azure",
                            linewidth=0.25,
                            edgecolor="xkcd:azure",
                        )
                    # plt.xlim(xmax=max_av)
                    # pdb.set_trace()

    plt.suptitle(r"Best-fit $A_V$ for each pixel", fontsize=40)

    pp.savefig()

    # -----------------
    # histograms of AVs with a chi2 cut
    # -----------------

    if len(chi2_plot) > 0:
        print("")
        print("Making Av Histograms with chi^2 cut")
        print("")

    for chi2_cut in chi2_plot:

        # set up figure
        fig = plt.figure(figsize=(x_dimen * 2, y_dimen * 2))

        for i in tqdm(range(y_dimen), desc="y pixels"):
            for j in tqdm(range(x_dimen), desc="x pixels"):
        # for i in [0]:
        #    for j in [12]:

                if nstars_image[i, j] > 20:

                    # set up the subplot
                    plt.subplot(y_dimen, x_dimen, (y_dimen - i - 1) * (x_dimen) + j + 1)

                    # make a histogram
                    if best_av[i][j] != []:
                        if not log_scale:
                            plot_av = np.array(best_av[i][j])[
                                np.array(best_av_chi2[i][j]) < chi2_cut
                            ]
                        if log_scale:
                            plot_av = np.log10(
                                np.array(best_av[i][j])[
                                    np.array(best_av_chi2[i][j]) < chi2_cut
                                ]
                            )
                        if len(plot_av) != 0:
                            h = plt.hist(
                                plot_av,
                                bins=bins.size,
                                range=(uniq_av[0] - gap / 2, uniq_av[-1] + gap / 2),
                                facecolor="xkcd:azure",
                                linewidth=0.25,
                                edgecolor="xkcd:azure",
                            )

        plt.suptitle(
            r"Best-fit $A_V$ for each pixel, but only using sources with $\chi^2 < $"
            + str(chi2_cut),
            fontsize=40,
        )

        pp.savefig()

    # close PDF figure
    pp.close()
