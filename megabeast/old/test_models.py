import pytest

from megabeast.mbsettings import mbsettings
from megabeast.singlepop_dust_model import MB_Model

fd_model = {
    "Av": {
        "name": "gaussian",
        "varnames": ["mean", "sigma"],
        "varinit": [1.0, 0.25],
        "prior": {
            "name": "flat",
            "var_minmax": [[0.005, 5.0], [0.05, 1.0]],
        },
    },
    "Rv": {
        "name": "gaussian",
        "varnames": ["mean", "sigma"],
        "varinit": [3.1, 0.25],
        "prior": {
            "name": "flat",
            "var_minmax": [[2.0, 6.0], [0.05, 1.0]],
        },
    }
}

models = [fd_model, fd_model]


@pytest.mark.parametrize("model", models)
def test_lnprior(model):
    """
    Test that the lnprior handles the defined prior types
    """
    priortypes = ["fixed", "flat"]

    # setup params
    params = mbsettings()
    params.fd_model = model
    mod = MB_Model(params)
    for cprior in priortypes:
        assert mod.lnprior(mod.start_params()) == 0.0, "test"
