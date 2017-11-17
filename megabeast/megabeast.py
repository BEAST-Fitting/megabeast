
# setup the megabeast model including defining the priors
#   - dust distribution model
#   - stellar populations model (later)

# use nstars image to setup for each pixel

# read in any info that is needed by all the pixels

# loop over the pixels with non-zero entries in the nstars image

# for each pixel
#    read in the sparse likelihoods and needed completeness
#    fit the megabeast model
#       - best fit via 'standard' fitter
#       - detailed likelhood via MCMC
#    save the fit info
#       - best fit
#       - megabeast parameter 1D pPDFs
#       - MCMC chain

# output the saved results as images and pixel specific files as needed
