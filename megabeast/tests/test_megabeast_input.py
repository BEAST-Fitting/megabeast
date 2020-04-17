import pkg_resources

from megabeast.read_megabeast_input import read_megabeast_input, required_params


def test_megabeast_input():
    """
    Test that the megabeast input parameter file reading works.
    """

    data_path = pkg_resources.resource_filename("megabeast", "examples/")
    tfilename = f"{data_path}/megabeast_input.txt"

    a = read_megabeast_input(tfilename)

    assert isinstance(a, dict), "result should be a dictionary"

    for tparam in required_params:
        assert tparam in a.keys(), "required parameter not present in file"
