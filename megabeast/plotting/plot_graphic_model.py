from graphviz import Digraph
import numpy as np

__all__ = ["create_graphic_model", "plot_graphic_model"]


def create_graphic_model(nodes, edges, gtype):
    """
    Create a graphic model given nodes and edges

    Parameters
    ----------
    nodes : dict
        for each node {key, text, math}
    edges : dict
        for each edge {key, text, math}
    gtype : str [default="text"]
        "text" for a verbose version, "math" for a compact version
    """
    mod = Digraph()
    if gtype == "math":
        tindx = 1
    else:
        tindx = 0
    for ckey in nodes.keys():
        if ckey == "Like":
            cstyle = "filled"
        else:
            cstyle = None
        mod.node(ckey, nodes[ckey][tindx], style=cstyle)

    for ckey in edges.keys():
        for cval in np.atleast_1d(edges[ckey]):
            mod.edge(ckey, cval)

    return mod


def plot_graphic_model(gtype="text", savefig="png"):
    """
    Plot the graphical model of the BEAST.

    Parameters
    ----------
    gtype : str [default="text"]
        "text" for a verbose version, "math" for a compact version
    savefig : str
        set to the file extension of desired file to save image of model
    """

    nodes = {
        "pIMF": ("IMF prior", "pIMF"),
        "IMF": ("Initial\nMass\nFunction", "<IMF(M<SUB>l</SUB>, M<SUB>h</SUB>, slope(s))>"),
        "pSFH": ("SFH prior", "pSFH"),
        "SFH": ("Star\nFormation\nHistory", "SFH(t)"),
        "pAMR": ("AMR prior", "pAMR"),
        "AMR": ("Age\nMetallicity\nRelation", "AMR(Z, t)"),
        "pDP": ("Distance prior", "pD"),
        "DP": ("Distance", "D(d)"),
        "pFC": ("prior\nForeground\nDust", "<pFD>"),
        "FC": ("Foreground\nDust", "<FD(A(V), R(V), f<SUB>A</SUB>)>"),
        "pGC": ("prior\nInternal\nDust", "<pID>"),
        "GC": ("Internal\nDust", "<ID(A(V), R(V), f<SUB>A</SUB>)>"),
        "M": ("mass\nM", "M"),
        "t": ("age\nT", "t"),
        "Z": ("metallicity\nZ", "Z"),
        "d": ("distance\nd", "d"),
        "Av": ("dust column\nA(V)", "A(V)"),
        "Rv": ("grain size\nR(V)", "R(V)"),
        "fA": ("<f<SUB>A</SUB>>", "<f<SUB>A</SUB>>"),
        "C": ("Completeness", "C(&theta;)"),
        "Cont": ("Contaminants\n&alpha;", "&alpha;"),
        "Phys": ("MegaBEAST\nPhysics+Observation Model", "p(&theta;|&phi;)"),
        "Like": ("BEAST\nLikelihoods", "BEAST\nL(&theta;)"),
    }

    edges = {
        "pIMF": "IMF",
        "IMF": "M",
        "pSFH": "SFH",
        "SFH": ("M", "t"),
        "AMR": ("Z", "t"),
        "pAMR": "AMR",
        "pDP": "DP",
        "DP": "d",
        "pFC": "FC",
        "FC": ("Av", "Rv", "fA"),
        "pGC": "GC",
        "GC": ("Av", "Rv", "fA"),
        "M": "Phys",
        "t": "Phys",
        "Z": "Phys",
        "d": "Phys",
        "Av": "Phys",
        "Rv": "Phys",
        "fA": "Phys",
        "C": "Phys",
        "Phys": "Like",
        "Cont": "Like",
    }

    beast = create_graphic_model(nodes, edges, gtype)

    beast.render(f"megabeast-graphic-{gtype}", format=savefig)


if __name__ == "__main__":
    plot_graphic_model("text")
    plot_graphic_model("math")
