from megabeast.mbsettings import mbsettings


def test_basic():
    mbparams = mbsettings()

    assert isinstance(mbparams, mbsettings)
