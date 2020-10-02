
##################
MegaBEAST Settings
##################

Physicsmodel
============

An example of the physicsmodel settings for a population of stars with one foreground dust cloud
and one internal dust cloud.  The internal cloud is assumed to be in an infinite
place in the center of populations of stars.  Thus, 1/2 the stars are in front and the other half
behind the internal cloud.  And all the stars are behind the foreground cloud.

The possible models are tabulated in the
`BEAST priors docs <https://beast.readthedocs.io/en/latest/beast_priors.html>`_.

::

  # stellar population model
  stellar_model = {
    "name": "bins_histo",
    "logages": [6.0, 7.0, 8.0, 9.0, 10.0],  # units are log(years)
    "imf_model": {"name": "kroupa"},
    "sfh_model": {
        "name": "flat",
        "init": [1e-5, 1e-5, 1e-5, 1e-5, 1e-5],  # units are M_sun/year
        "prior": {
            "name": "flat",
            "min": [0.0, 0.0, 0.0, 0.0, 0.0],
            "max": [1e-4, 1e-4, 1e-4, 1e-4, 1e-4],
        },
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
    "Av_model": {
        "name": "log-normal",
        "m_init": 0.2,
        "s_init": 0.1,
        "prior": {
            "name": "flat",
            "m_min": 0.005,
            "m_max": 0.5,
            "s_min": 0.05,
            "s_max": 0.2,
        },
    },
    "Rv_model": {"name": "flat", "m_init": 3.1, "prior": {"name": "fixed"}},
    "fA_model": {"name": "flat", "m_init": 1.0, "prior": {"name": "fixed"}},
  }

  # internal dust cloud
  id_model = {
    "Av_model": {
        "name": "log-normal",
        "m_init": 1.0,
        "s_init": 0.5,
        "prior": {
            "name": "flat",
            "m_min": 0.1,
            "m_max": 10.0,
            "s_min": 0.1,
            "s_max": 1.0,
        },
    },
    "Rv_model": {
        "name": "log-normal",
        "m_init": 1.0,
        "s_init": 0.5,
        "prior": {
            "name": "flat",
            "m_min": 0.1,
            "m_max": 10.0,
            "s_min": 0.1,
            "s_max": 1.0,
        },
    },
    "fA_model": {
        "name": "log-normal",
        "m_init": 1.0,
        "s_init": 0.5,
        "prior": {
            "name": "flat",
            "m_min": 0.1,
            "m_max": 10.0,
            "s_min": 0.1,
            "s_max": 1.0,
        },
    },
  }
