#############################
Simulations for the MegaBEAST
#############################

**********
Background
**********

To test that the MegaBEAST is behaving as expected, it is important to simulate data and evaluate if the results match the input.  Currently, code is implemented to create a lognormal distribution of A_V values that can be run through the MegaBEAST.  More sophisticated inputs, as well as additional types of simulations, are still necessary.

The code `simulate_av.py` will create all files necessary to run the MegaBEAST.  The inputs are:

  * `beast_seds_filename`: Name of the path+file that contains the BEAST SED grid
  * `av_lognorm`: parameters for the underlying lognormal A_V distribution
  * `output_label`: a string used for making output folders and file names
  * `image_dimen` (optional): dimensions of the fake sky image (default is 10x10)
  * `nstar_per_pix` (optional): number of stars in each pixel (default is 50)

The code will create an `nstars` image with the specified dimensions.  Each pixel will be populated with a number of stars drawn from a poisson distribution with mean of `nstar_per_pix`.
In each pixel, each star will be assigned an A_V drawn from the input lognormal distribution. The error on A_V is currently assumed to be 0.2 (that will be made more generalizable in the future).
The log-likelihoods (`lnp`) are computed at each A_V in the BEAST grid.  Each star is assumed to have a completeness of 1 (which will also be more generalized in the future).

The code creates a folder with the spatially organized log likelihoods, an `nstars` image, and a placeholder noisemodel.  It also writes out a MegaBEAST input file, which can be directly input into the MegaBEAST.

*******
Example
*******

This example creates an 8x10 image in which each pixel has 50 stars.  The assumed underlying lognormal distribution of A_V has a peak at A_V=1 and sigma=0.5.  The BEAST grid has A_V spacing of 0.01.

.. code-block:: shell

  >>> from simulate_av import simulate_av
  >>> from megabeast import megabeast
  >>> from plot_input_data import plot_input_data
  >>> simulate_av('someproject_beast_seds.grid.hd5',
                  {'max_pos':1, 'sigma':0.5, 'N':500},
		  'mb-simulation',
		  image_dimen=[8,10],
		  nstar_per_pix=50)
  >>> megabeast('megabeast_input_mb-simulation.txt')
  >>> plot_input_data('megabeast_input_mb-simulation.txt', log_scale=True)

This creates the simulated files, runs the MegaBEAST, and makes plots.

