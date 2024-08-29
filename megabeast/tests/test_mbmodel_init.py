import astropy.units as u
from megabeast.mbmodel import MBModel


def test_mbmodel_init():

    # stellar population model
    stellar_model = {
        "logA": {  # star formation history SFH
            "name": "bins_histo",
            "x": [6.0, 7.0, 8.0, 9.0, 10.0],  # units are log(years)
            "varnames": ["values"],
            "varinit": [[1e-8, 1e-8, 1e-8, 1e-8, 1e-8]],  # units are M_sun/year
            "prior": {
                "name": "flat",
                "minmax": [[0.0, 0.0, 0.0, 0.0, 0.0], [1e-3, 1e-3, 1e-3, 1e-3, 1e-3]],
            },
        },
        "M_ini": {  # initial mass function
            "name": "salpeter",
            "varnames": ["slope"],
            "varinit": [2.35],
            "prior": {
                "name": "flat",
                "minmax": [[2.0, 3.0]],
            },
        },
        "Z": {
            "name": "flat",
        },
        "distance": {
            "name": "absexponential",
            "varnames": ["dist0", "tau", "amp"],
            "varinit": [60.0 * u.kpc, 5.0 * u.kpc, 1.0],
            "prior": {
                "name": "flat",
                "minmax": [[50.0, 3.0, 0.9], [70.0, 7.0, 1.1]],
            },
        },
    }

    # foreground dust cloud
    dust_model = {
        "Av": {
            "name": "lognormal",
            "varnames": ["mean", "sigma"],
            "varinit": [1.0, 0.1],
            "prior": {
                "name": "flat",
                "minmax": [[0.005, 5.0], [0.05, 1.0]],
            },
        },
        "Rv": {
            "name": "lognormal",
            "varnames": ["mean", "sigma"],
            "varinit": [3.1, 0.25],
            "prior": {
                # "name": "fixed",
                "name": "flat",
                "minmax": [[2.0, 6.0], [0.05, 1.0]],
            },
        },
        "f_A": {
            "name": "lognormal",
            "varnames": ["mean", "sigma"],
            "varinit": [1.0, 0.25],
            "prior": {
                "name": "flat",
                "minmax": [[0.0, 1.0], [0.05, 0.5]],
            },
        },
    }

    mod = MBModel(stellar_model, dust_model)

    assert isinstance(mod.star_model, dict)
    assert isinstance(mod.dust_model, dict)

    assert isinstance(mod.physics_model, dict)

    assert isinstance(mod.compute_N_stars, bool)
    assert isinstance(mod.compute_massmult, bool)

    # parameters expected
    expparam = ["logA", "M_ini", "Z", "distance", "Av", "Rv", "f_A"]
    for cparam in mod.physics_model.keys():
        assert cparam in expparam
        for csubparam in ["name", "model"]:
            assert (
                csubparam in mod.physics_model[cparam].keys()
            ), f"required {csubparam} not found in {cparam}"
