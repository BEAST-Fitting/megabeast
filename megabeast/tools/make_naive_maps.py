import math

import argparse
import numpy as np
import h5py
import itertools as it
from astropy.io import fits
from astropy.table import Table, vstack

from beast.tools.create_background_density_map import (
    calc_nx_ny_from_pixsize,
    make_wcs_for_map,
    get_pix_coords,
)


__all__ = ["create_naive_maps"]


def create_naive_maps(stats_filename,
                      pix_size=10.0,
                      verbose=False,
                      median=False,
                      chi2mincut=None,
                      weigh_by_av=False):
    """
    Make the naive maps by directly averaging the BEAST results for all the
    stars in each pixel.  Does not account for completeness, hence naive maps!

    Parameters
    ----------
    stats_filename : string or list of strings
       name(s) of the catalog(s) of BEAST results

    pix_size : float (default=10)
       size of pixels/regions in arcsec

    median : bool (default=False)
       calculate the median of the BEAST results (instead of the mean)

    chi2mincut : int (default=None)
       place a chi2min cut on the BEAST results. Gives the max threshold.

    weigh_by_av : bool (default=False)
        weigh R(V) and f_A by A(V) to determinе R(V) and f_A of the total column
        of dust in a pixel (as opposed to finding a simple average across a pixel)
    """

    # type of statistic (make a commandline parameter later)
    #   remember to add to output filenames
    stat_type = "Exp"

    # read in the full brick catalog
    if isinstance(stats_filename, str):
        stats_filename = [stats_filename]
    cat = Table.read(stats_filename[0])
    if len(stats_filename) > 1:
        for fname in stats_filename[1:]:
            tcat = Table.read(fname)
            cat = vstack([cat, tcat])

    if chi2mincut:
        cat = cat[cat['chi2min'] < chi2mincut]

    # make RA/Dec grid
    ra = cat["RA"]
    dec = cat["DEC"]
    pixsize_degrees = pix_size / 3600
    n_x, n_y, ra_delt, dec_delt = calc_nx_ny_from_pixsize(cat, pixsize_degrees)
    # the ra spacing needs to be larger, as 1 degree of RA ==
    # cos(DEC) degrees on the great circle
    ra_grid = ra.min() + ra_delt * np.arange(0, n_x + 1, dtype=float)
    dec_grid = dec.min() + dec_delt * np.arange(0, n_y + 1, dtype=float)

    # generate the wcs info for the output FITS files
    w = make_wcs_for_map(ra_grid, dec_grid)

    # get the pixel coordinates for each source
    pix_x, pix_y = get_pix_coords(cat, w)
    # import pdb; pdb.set_trace()

    # for ease of checking the bin, set x/y coords to integers
    x = np.floor(pix_x)
    y = np.floor(pix_y)

    # setup arrary to store summary stats per pixel
    sum_stats = ["Av", "Rv", "f_A", "logT", "M_act", "logA"]
    print("sum_stats", sum_stats)
    n_sum = len(sum_stats)
    summary_stats = np.zeros((n_y + 1, n_x + 1, n_sum + 1), dtype=float)
    summary_sigmas = np.zeros((n_y + 1, n_x + 1, n_sum), dtype=float)
    values_foreach_pixel = {
        cur_stat: {(i, j): [] for i in range(n_x + 1) for j in range(n_y + 1)}
        for cur_stat in sum_stats
    }

    # loop through the pixels and generate the summary stats
    for i in range(n_x + 1):
        for j in range(n_y + 1):
            (tindxs,) = np.where((x == i) & (y == j))
            if len(tindxs) > 0:
                summary_stats[j, i, n_sum] = len(tindxs)
                if verbose:
                    print(i, j, len(tindxs))
                for k, cur_stat in enumerate(sum_stats):
                    values = cat[cur_stat + "_" + stat_type][tindxs]

                    # weigh R(V) and f_A by A(V)
                    if weigh_by_av and ("Rv" in cur_stat or "f_A" in cur_stat):
                        # get Av values
                        av_values = cat["Av" + "_" + stat_type][tindxs]
                        values_foreach_pixel[cur_stat][i, j] = values / av_values
                        weights = av_values
                    else:
                        values_foreach_pixel[cur_stat][i, j] = values
                        weights = None

                    if median:
                        summary_stats[j, i, k] = np.median(values)
                        xabs = abs(values - np.median(values)) ** 2.
                        summary_sigmas[j, i, k] = np.sqrt(np.median(xabs)) / math.sqrt(len(values))
                    else:
                        summary_stats[j, i, k] = np.average(values, weights=weights)
                        summary_sigmas[j, i, k] = np.std(values, ddof=1) / math.sqrt(len(values))

    master_header = w.to_header()
    # Now, write the maps to disk
    for k, cur_stat in enumerate(sum_stats):
        map_name = stats_filename[0].replace("stats.fits", "map_" + cur_stat + "_" + str(pix_size) + "arcsec.fits")
        hdu = fits.PrimaryHDU(summary_stats[:, :, k], header=master_header)
        hdu.writeto(map_name, overwrite=True)

        sigma_name = map_name.replace("map", "map_sigma")
        hdu_sigma = fits.PrimaryHDU(summary_sigmas[:, :, k], header=master_header)
        hdu_sigma.writeto(sigma_name, overwrite=True)

    hdu = fits.PrimaryHDU(summary_stats[:, :, n_sum], header=master_header)
    hdu.writeto(stats_filename[0].replace("stats.fits", "npts.fits"), overwrite=True)

    # And store all the values in HDF5 format
    values_name = stats_filename[0].replace("stats.fits", "values_per_pixel.hd5")
    f = h5py.File(values_name, "w")
    dt = h5py.special_dtype(vlen=float)
    for cur_stat in sum_stats:
        dset = f.create_dataset(cur_stat, (n_x, n_y), dtype=dt)
        for i, j in it.product(range(n_x), range(n_y)):
            dset[i, j] = values_foreach_pixel[cur_stat][i, j]


if __name__ == "__main__":

    # command line params to specify the run directory
    #   and any other needed parameters

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "stats_filename",
        metavar="fname",
        type=str,
        nargs="+",
        help="Filename(s) of the stats",
    )
    parser.add_argument(
        "--pix_size", default=10.0, type=float, help="pixel scale [arcsec]"
    )
    parser.add_argument(
        "--verbose", default=False, type=bool, help="print pixel indices as a check"
    )
    parser.add_argument(
        "--median", default=False, type=bool, help="find the median of the values"
    )
    parser.add_argument(
        "--chi2mincut", default=False, type=int, help="max chi2min threshold to place on results"
    )
    parser.add_argument(
        "--weigh_by_av", default=False, type=bool, help="weigh R(V) or f_A by A(V)"
    )
    args = parser.parse_args()

    # call the function
    create_naive_maps(args.stats_filename, pix_size=args.pix_size)
