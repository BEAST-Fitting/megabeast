####################
MegaBEAST Motivation
####################

The overarching motivation for the MegaBEAST is to fit for the ensemble
stellar and dust extinction parameters of observations of a set of stars
in a galaxy.
One key part of the MegaBEAST ensemble model is to include the effects
of completeness from the observations where these effects are calculated from
Artificial Star Test (ASTs).
Example ensemble stellar parameters are initial mass function (IMF) and
star formation history (SFH).
Example ensemble dust extinction parameters are maps of dust column A(V) and
average grain size R(V).
The MegaBEAST is a hierarchical Bayesian model that uses the results of the
BEAST Bayesian model.

The BEAST provides fits to photometric SEDs of individual stars
giving stellar and dust extinction parameters.
The BEAST results include various statistics for the BEAST primary and
derived parameters from best fits to the full n-dimensional posterior.
The BEAST does not account for completeness as it only fits individual stars.

The MegaBEAST takes advantage of all the BEAST work by using the
BEAST results for each star (specifically the n-dimensional posterior).
Thus, the MegaBEAST needs to be able to predict the stellar and dust
properties of ensembles of extinguished stars in the BEAST parameters only.
The MegaBEAST does need to include the observational effect of a finite
depth survey of stars, hence the ensemble model includes a model of the
completeness.

Goals
=====

The goal of the MegaBEAST is to derive maps of stellar and dust extinction
parameters from multi-band resolved star surveys of galaxies.

The mapped ensemble parameters include (but are not limited to):

- star formation history
- initial mass function
- mass-metallicity relationship
- column A(V)
- average grain size R(V)
- grain composition measure f_A
