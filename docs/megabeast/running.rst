#####################
Running the MegaBEAST
#####################

The MegaBEAST is specifically setup to use the output of BEAST
runs.  Tt uses the results ordered spatially to have
all the stars in a region split in to smaller pixels.  For example,
the BEAST results for a PHAT brick will have many files with each
file giving the results for a 5x5 arcsec^2 pixel.  The MegaBEAST uses
the BEAST nD pPDFs, usually in the form of a small subset of the
full pPDFs points sampled from the region with significant probability.

The files the MegaBEAST expects in the directory have filenames in the
format of projectname_XX_XX_lnp.hd5 where the XX and YY values are the
pixel coordinates x and y.  The translation between the x,y values and
ra,dec are given in the projectname_nstars.fits file that was created
when the BEAST spatial reordering was done.

More information on the BEAST spatial reordering can be found in the
`BEAST documentation <http://beast.readthedocs.io/en/latest/workflow.html#post-processing>`_.

Once the spatially reordered BEAST data is ready, the MegaBEAST is run
from the commandline.

.. code-block:: shell

  $ python megabeast.py projectname

The options can be found with the usual command.

.. code-block:: shell

  $ python megabeast.py --help
  usage: megabeast.py [-h] [--min_for_fit MIN_FOR_FIT] [-v] projectname

  positional arguments:
    projectname           project name to use (basename for files)

  optional arguments:
    -h, --help            show this help message and exit
    --min_for_fit MIN_FOR_FIT
                          minimum number of stars need in a pixel for fit
    -v, --verbose         verbose output

The results of the MegaBEAST are the best fit values for the ensemble model
for the A(V) distribution.  Currently, only the A(V) ensemble model is
implemented.  The output filenames have the format
projectname_param_bestfit.fits where the allowed values of param are
N1, N2, AV1, AV2, sigma1, and sigma2.  There are two lognormal components in
the A(V) ensemble model where the lognormal paramters are N (number of stars),
AV (peak A(V)), and sigma (width of lognormal).
