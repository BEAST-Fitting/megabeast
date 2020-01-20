from astropy.table import Table
import ast


def read_megabeast_input(input_file):
    """
    Read in the file with MegaBEAST settings
    """

    # read everything in
    input_data = Table.read(
        input_file, format="ascii.no_header", delimiter="=", comment="#"
    )

    # parse it
    input_vars = {}

    for i in range(len(input_data)):
        try:
            input_vars[input_data["col1"][i]] = ast.literal_eval(input_data["col2"][i])
        except Exception:
            input_vars[input_data["col1"][i]] = str(input_data["col2"][i])

    # verify that the necessary parameters are present
    required_params = [
        "beast_seds_filename",
        "beast_noise_filename",
        "av_prior_model",
        "rv_prior_model",
        "fA_prior_model",
        "nstars_filename",
        "lnp_file_prefix",
        "projectname",
        "fit_param_names",
        "min_for_fit",
    ]
    for param in required_params:
        if param not in input_vars.keys():
            raise ValueError("MegaBEAST input file is missing a parameter: " + param)

    # return it
    return input_vars
