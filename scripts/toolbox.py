"""
Toolbox of functions for bayes_mod

author: Andrew McCluskey (andrew.mccluskey@diamond.ac.uk)
"""

import periodictable as pt
import numpy as np
import scipy.constants as const
from refnx.analysis import Transform, Objective, GlobalObjective
from refnx.reflect import ReflectModel, SLD


SCALES = [0.19, 0.22, 0.21, 0.15, 0.15, 0.18, 0.15]


def get_scattering_length(component, xray=False):
    """
    Determine the scattering length for a dictionary of atoms

    Parameters
    ----------
    component : dict
        Dictionary containing element (isotope) types and numbers of each
    xray : bool
        If X-ray scattering factor to be evaluated

    Returns
    -------
    scattering_length : complex
        Scattering length for the given component
    """
    if xray:
        scattering_length = 0 + 0j
        cre = const.physical_constants["classical electron radius"][0] * 1e10
        for key in component:
            scattering_length += np.multiply(
                pt.elements.symbol(key).xray.scattering_factors(energy=12)[0],
                cre * component[key],
            )
            scattering_length += (
                np.multiply(
                    pt.elements.symbol(key).xray.scattering_factors(energy=12)[
                        1
                    ],
                    cre * component[key],
                )
                * 1j
            )
    else:
        scattering_length = 0 + 0j
        for key in component:
            scattering_length += (
                pt.elements.symbol(key).neutron.b_c * component[key]
            )
            if pt.elements.symbol(key).neutron.b_c_i:
                inc = pt.elements.symbol(key).neutron.b_c_i
            else:
                inc = 0
            scattering_length += inc * 1j * component[key]
        scattering_length *= 1e-5
    return scattering_length


def get_objectives(num_contrasts, models, datasets):
    """
    Create global objective

    Parameters
    ----------
    num_contrasts : int
        Number of contrasts
    models : list of refnx.reflect.ReflectModel
        ReflectModels list
    datasets : refnx.dataset.ReflectDataset
        Experimental datasets list

    Returns
    -------
    g_objective : refnx.analysis.GlobalObjective
        Global objective for analysis
    """
    objectives = []
    for i in range(num_contrasts):
        objectives.append(
            Objective(models[i], datasets[i], transform=Transform("YX4"))
        )
    return GlobalObjective(objectives)


def get_models(num_contrasts, structures, datasets):
    """
    Create reflect models

    Parameters
    ----------
    num_contrasts : int
        Number of contrasts
    structures : list of refnx.reflect.Structure
        Monolayer structure list
    datasets : refnx.dataset.ReflectDataset
        Experimental datasets list

    Returns
    -------
    models : list of refnx.reflect.ReflectModel
        ReflectModels list
    """
    models = []
    for i in range(num_contrasts):
        models.append(ReflectModel(structures[i]))
        models[i].scale.setp(SCALES[i], vary=False)
        models[i].bkg.setp(datasets[i].y[-1], vary=False)
    return models


def get_structures(contrasts, lipids):
    """
    Create monolayer structures

    Parameters
    ----------
    contrasts : list
        List of contrasts
    lipids : refnx compatibles class
        Lipids objects

    Returns
    -------
    structures : list of refnx.reflect.Structure
        Monolayer structure list
    """
    air = SLD(0, "air")
    structures = []
    for i, cont in enumerate(contrasts):
        if cont[-3:] == "d2o":
            water = SLD(6.35, "d2o")
        elif cont[-4:] == "acmw":
            water = SLD(0.0, "acmw")
        structures.append(air(0, 0) | lipids[i] | water(0, 3.3))
    return structures


def get_lipids(contrasts, b_val, model):
    """
    Create lipid objects

    Parameters
    ----------
    contrasts : list
        List of contrasts
    b_val : list
        List of scattering lengths for lipid components
    model : refnx compatibles class
        Particular class to be used for lipids

    Returns
    -------
    lipids : list of refnx compatibles class
        Lipids object list
    """
    lipids = []
    for i, cont in enumerate(contrasts):
        lipids.append(model(b_val[i], name=cont))
    return lipids
