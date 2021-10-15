.. _mb-settings:

##################
MegaBEAST Settings
##################

.. note::
   The format of the settings file is still evolving.

BEAST Info
==========

The base name of the files is needed as this defines the BEAST physicsmodel
(basename_sed.grid.hd5), observationmodel (basname__noisemodel.grid.hd5),
and likelihoods (basename_lnp.hd5) files that are  used for the fitting.

::

   # BEAST results files
   beast_base = "beast_sim/scylla_sim/scylla_sim"

Physicsmodel
============

An example of the physicsmodel settings for a population of stars with one foreground dust cloud.
The internal cloud is assumed to be in an infinite
All the stars are behind the foreground cloud.

The possible models are tabulated in the
`BEAST priors docs <https://beast.readthedocs.io/en/latest/beast_priors.html>`_.

::

  # stellar population model
  stellar_model = {
      "name": "bins_histo",
      "x": [6.0, 7.0, 8.0, 9.0, 10.0],  # units are log(years)
      "logA": {  # star formation history SFH
          "name": "bins_histo",
          "varnames": ["values"],
          "varinit": [[1e-8, 1e-8, 1e-8, 1e-8, 1e-8]],  # units are M_sun/year
          "prior": {
              "name": "flat",
              "minmax": [[0.0, 0.0, 0.0, 0.0, 0.0], [1e-3, 1e-3, 1e-3, 1e-3, 1e-3]],
          },
      },
      "M_ini": {  # initial mass function
          "name": "kroupa",
          "varnames": ["slope"],
          "varinit": [2.35],
          "prior": {
            "name": "fixed",
            "minmax": [[2.0, 3.0]],
          }
      },
      "amr_model": {
          "name": "bins_histo",
          "init": [-0.3, -0.3, -0.3, -0.3, -0.3],
          "prior": {
              "name": "flat",
              "min": [-2.1, -2.1, -2.1, -2.1, -2.1],
              "max": [0.0, 0.0, 0.0, 0.0, 0.0],
          },
      },
      "dist_model": {
          "name": "Gaussian",
          "m_unit": "kpc",
          "m_init": 60.0,  # units are
          "s_init": 5.0,
          "d_prior": {"name": "flat", "min": 50.0, "max": 70.0},
          "s_prior": {"name": "flat", "min": 1.0, "max": 10.0},
      },
  }

  # foreground dust cloud
  fd_model = {
      "Av": {
          "name": "lognormal",
          "varnames": ["mean", "sigma"],
          "varinit": [0.5, 0.25],
          "prior": {
              "name": "flat",
              "minmax": [[0.005, 5.0], [0.05, 1.0]],
          },
      },
      "Rv": {
          "name": "lognormal",
          "varnames": ["mean", "sigma"],
          "varinit": [4.0, 0.25],
          "prior": {
              "name": "flat",
              "minmax": [[2.0, 6.0], [0.05, 1.0]],
          },
      },
      "fA": {
          "name": "lognormal",
          "varnames": ["mean", "sigma"],
          "varinit": [1.0, 0.25],
          "prior": {
              "name": "fixed",
              "minmax": [[0.0, 1.0], [0.05, 0.5]],
          },
      }
  }
