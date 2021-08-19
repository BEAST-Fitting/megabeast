import numpy as np

from beast.tools.read_beast_data import read_lnp_data

__all__ = ["get_likelihoods"]


def get_likelihoods(ppdf_file, beast_model_data):
    """
    Read in the saved BEAST sparse posterior PDFs and divide by the BEAST
    priors to get the spare liklilhoods.

    Parameters
    ----------
    ppdf_file : string
        filename of the saved BEAST posterior PDFs

    Returns
    -------
    lnp_data : dictonary
       contains arrays of the likelihood values (vals) and indices to the BEAST model grid (indxs)
    """
    # BEAST saves posterior PDFs labeled as log(pPDF)
    lnpdata = read_lnp_data(ppdf_file)

    # divide by the BEAST prior weights to return to recover the likelihoods
    n_lnps, n_stars = lnpdata["indxs"].shape
    for i in range(n_stars):
        indxs = lnpdata["indxs"][:, i]
        gmask = np.isfinite(lnpdata["vals"][:, i])
        lnpdata["vals"][gmask, i] = (
            np.exp(lnpdata["vals"][gmask, i])
            / beast_model_data["prior_weight"][indxs[gmask]]
        )

    return lnpdata
