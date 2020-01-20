# system
import os
import glob

# other packages
import numpy as np
import scipy.stats
from astropy.io import fits
import h5py
from astropy import wcs

# beast
from beast.physicsmodel.prior_weights_dust import _lognorm
from beast.fitting.fit import save_lnp


def simulate_av(
    beast_seds_filename,
    av_lognorm,
    output_label,
    av_lognorm2=None,
    image_dimen=[10, 10],
    nstar_per_pix=50,
):
    """
    Create a set of fake data in which A_V follows a log-normal distribution, which can then be run through the megabeast

    Parameters
    ----------
    beast_seds_filename : string
        Name of the file that contains the BEAST SED grid

    av_lognorm : dict
        Dictionary with parameters for lognormal distribution (keys = 'max_pos', 'sigma', 'N')

    output_label : string
        label to use for folders and file names

    av_lognorm2 : dict (default=None)
        if set (identical formatting to av_lognorm), makes A_V a double lognormal distribution

    image_dimen : list of ints (default=[10,10])
        dimensions (nx by ny) of the fake image

    nstar_per_pix : int (default=50)
        number of stars per pixel (this will vary between pixels using a Poisson distribution)

    """

    # check that the directory exists
    if not os.path.exists("./" + output_label + "/"):
        os.makedirs("./" + output_label + "/")

    # read in A_V info from the SED grid data
    beast_seds_hdf = h5py.File(beast_seds_filename, "r")
    beast_av_list = np.unique(beast_seds_hdf["grid"]["Av"])
    # get first index of each A_V
    sed_grid_length = len(beast_seds_hdf["grid"]["Av"])
    beast_av_ind = [
        int(sed_grid_length / len(beast_av_list) * i) for i in range(len(beast_av_list))
    ]

    # create an nstars image
    nstars_hdu = setup_nstars_image(image_dimen, nstar_per_pix, output_label)

    # create lnp files for each pixel
    # - first, delete any existing ones (otherwise new stuff is appended)
    for f in glob.glob("./" + output_label + "/" + output_label + "_*_*_lnp.hd5"):
        os.remove(f)
    # - now make new ones
    setup_lnp_files(
        nstars_hdu,
        av_lognorm,
        beast_av_list,
        np.array(beast_av_ind),
        output_label,
        av_lognorm2=av_lognorm2,
    )

    # create a noise file (currently all stars get completeness=1)
    setup_noise_file(sed_grid_length, output_label)

    # create a megabeast input file
    setup_mb_input(beast_seds_filename, output_label)


def setup_mb_input(beast_seds_filename, output_label):
    """
    Write out a megabeast input file

    output_label is enough to construct all the necessary file names
    """

    file_lines = [
        "# -----------",
        "# BEAST files",
        "# -----------",
        "",
        "# files that contain the SEDs and noise model from the BEAST",
        "beast_seds_filename = " + beast_seds_filename,
        "beast_noise_filename = ./"
        + output_label
        + "/"
        + output_label
        + "_noisemodel.hd5",
        "",
        "# prior models for those",
        "av_prior_model = {'name': 'flat'}",
        "rv_prior_model = {'name': 'flat'}",
        "fA_prior_model = {'name': 'flat'}",
        "",
        "# fits file that contains the number of stars in each spatially-reordered pixel",
        "nstars_filename = ./" + output_label + "/" + output_label + "_nstars.fits",
        "",
        "# prefix for log likelihoods",
        "# file will be constructed as lnp_file_prefix_#_#_lnp.hd5, where #s are the spatially reordered pixel numbers",
        "lnp_file_prefix = ./" + output_label + "/" + output_label,
        "",
        "# ------------------",
        "# MegaBEAST settings",
        "# ------------------",
        "",
        "# project name to use for megabeast output files",
        "# megabeast outputs will be saved in projectname_megabeast folder",
        "projectname = " + output_label,
        "",
        "# parameters for the megabeast to fit",
        "fit_param_names = ['Av1', 'Av2', 'sigma1', 'sigma2', 'N12_ratio']",
        "",
        "# minimum number of stars need in a pixel for fit",
        "min_for_fit = 20",
        "",
    ]

    with open("megabeast_input_" + output_label + ".txt", "w") as mb_file:
        for line in file_lines:
            mb_file.write(line + "\n")


def setup_noise_file(sed_grid_length, output_label):
    """
    Create a noise file

    Currently, the megabeast only uses the completeness column, but the file
    also has bias and error columns.  For now, they'll all be set as:
    bias = 0
    completeness = 1
    error = 0

    Each array has shape (n_stars, n_filters) -> use n_filters = 1 for
    simplicity, since megabeast uses the max completeness
    """

    # create the file
    with h5py.File(
        "./" + output_label + "/" + output_label + "_noisemodel.hd5", "w"
    ) as f:

        # create the data sets
        f.create_dataset("bias", data=np.zeros((sed_grid_length, 1)))
        f.create_dataset("completeness", data=np.ones((sed_grid_length, 1)))
        f.create_dataset("error", data=np.zeros((sed_grid_length, 1)))


def setup_lnp_files(
    nstars_hdu, av_lognorm, av_gridpoints, av_ind, output_label, av_lognorm2=None
):
    """
    Create sparsely sampled lnp files for each pixel

    Parameters
    ----------
    nstars_hdu : hdu object
        hdu with the nstars image

    av_lognorm : dict
        Dictionary with parameters for lognormal distribution (keys = 'max_pos', 'sigma', 'N')

    av_gridpoints : array
        array of all the unique A_V values in the BEAST grid

    av_ind : list of ints
        indices in the BEAST grid for each A_V in av_gridpoints

    output_label : string
        label to use for folders and file names

    av_lognorm2 : dict (default=None)
        if set (identical formatting to av_lognorm), makes A_V a double lognormal distribution

    """

    # dimensions of image
    nx, ny = nstars_hdu.data.shape

    # create a lognormal distribution
    av_dist = np.linspace(av_gridpoints[0], av_gridpoints[-1], 500)
    temp = _lognorm(
        av_dist, av_lognorm["max_pos"], sigma=av_lognorm["sigma"], N=av_lognorm["N"]
    )
    if av_lognorm2 is not None:
        temp2 = _lognorm(
            av_dist,
            av_lognorm2["max_pos"],
            sigma=av_lognorm2["sigma"],
            N=av_lognorm2["N"],
        )
        temp += temp2
    lognorm_dist = temp / np.sum(temp)
    lognorm_cdf = np.cumsum(lognorm_dist)

    # iterate through pixels
    for i in range(nx):
        for j in range(ny):

            # number of stars in this pixel
            nstar = nstars_hdu.data[i, j]

            # initialize a list to save info for lnp file
            lnp_save_list = []

            # iterate through the stars in this pixel
            for n in range(nstar):

                # get an A_V from the lognormal distribution
                av = np.interp(np.random.random(1), lognorm_cdf, av_dist)

                # calculate probabilities
                # include errors -> sample from gaussian at each BEAST A_V grid point
                prob_list = scipy.stats.norm(av, 0.2).pdf(av_gridpoints)
                prob_list = prob_list / np.sum(prob_list)

                # make a list of the required quantities in the lnp file
                star_num = n
                idx = av_ind
                lnp = np.log(prob_list)
                chi2 = lnp / (-0.5)
                sed = [99]
                lnp_save_list.append([star_num, idx, lnp, chi2, sed])

            # save the lnp data
            # print('i='+str(i)+' j='+str(j)+' nstar='+str(nstar)+' len(lnp_list)='+str(len(lnp_save_list)))
            save_lnp(
                output_label + "/" + output_label + "_%i_%i_lnp.hd5" % (j, i),
                lnp_save_list,
                False,
            )


def setup_nstars_image(image_dimen, nstar_per_pix, output_label):
    """
    Create nstars image

    Returns an HDU with the image
    """

    # create WCS (exact values don't matter, but the megabeast code requires WCS)
    # - arbitrary settings
    dec_delt = 10.0 / 3600.0  # (10 arcsec)
    ra_delt = dec_delt
    min_ra, max_ra = 0, image_dimen[0] * ra_delt
    min_dec, max_dec = 0, image_dimen[1] * dec_delt
    # - WCS thingy
    w = wcs.WCS(naxis=2)
    w.wcs.crpix = np.asarray(image_dimen, dtype=float) / 2.0
    w.wcs.cdelt = [ra_delt, dec_delt]
    w.wcs.crval = np.asarray([(min_ra + max_ra), (min_dec + max_dec)]) / 2.0
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    # - turn it into a fits header
    master_header = w.to_header()

    # create image
    nstar_array = np.random.poisson(nstar_per_pix, tuple(image_dimen))

    # create HDU
    hdu = fits.PrimaryHDU(data=nstar_array, header=master_header)

    hdu.writeto(output_label + "/" + output_label + "_nstars.fits", overwrite=True)

    return hdu
