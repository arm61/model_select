"""
Calculate the reflectometry with respect to a set of contrast from the
two_layer model.

Copyright (c) Andrew R. McCluskey

Distributed under the terms of the MIT License

author: Andrew R. McCluskey (andrew.mccluskey@diamond.ac.uk)
"""

# pylint: disable=R0912

from refnx.dataset import ReflectDataset
import toolbox
import two_layer


DATA_DIR = "data"


def refl(values, contrasts, variables, logl=True):
    """
    Evalutate the reflectometry from two-layer phospholipid monolayer model.

    Args:
        values (array_like): An array of values for the parameters that are
            being varied.
        contrasts (list): A list of strings describing the contrasts being
            investigated.
        variables (dict): A dict describing the parameters to be varied and
            the index of this parameter in the values array.
        logl (bool, optional): Return just the logl value (if `True`) or the
            whole refnx.analysis.GlobalObjective object (if `False`).

    Returns:
        (float or refnx.analysis.GlobalObjective) Depending on the value of
            `logl`, either the logl of the input or the whole
            GlobalObjective is returned.
    """
    datasets = []
    num_contrasts = len(contrasts)
    for i in range(num_contrasts):
        datasets.append(
            ReflectDataset("{}/{}30.dat".format(DATA_DIR, contrasts[i]))
        )

    b_lipid = two_layer.get_b_lipid(contrasts)

    lipids = toolbox.get_lipids(contrasts, b_lipid, two_layer.TwoLayer)

    structures = toolbox.get_structures(contrasts, lipids)

    for i, contrast in enumerate(lipids):
        if "th" in variables:
            contrast.thick_h.setp(values[variables["th"]], vary=False)
        else:
            contrast.thick_h.setp(10.0, vary=False)
        if "mvh" in variables:
            contrast.mol_vol_h.setp(values[variables["mvh"]], vary=False)
        else:
            contrast.mol_vol_h.setp(339.5, vary=False)
        if "tt" in variables:
            contrast.thick_t.setp(values[variables["tt"]], vary=False)
        else:
            contrast.thick_t.setp(21, vary=False)
        if "phit" in variables:
            contrast.phi_t.setp(values[variables["phit"]], vary=False)
        else:
            contrast.phi_t.setp(1.0, vary=False)
        if "mvt" in variables:
            contrast.mol_vol_t.setp(values[variables["mvt"]], vary=False)
        else:
            contrast.mol_vol_t.setp(850.4, vary=False)
        contrast.phi_h.constraint = (
            contrast.thick_t
            * contrast.mol_vol_h * contrast.phi_t
            / (contrast.mol_vol_t * contrast.thick_h)
        )
        if "rough" in variables:
            structures[i][-1].rough.setp(
                values[variables["rough"]],
                vary=False,
            )
        else:
            structures[i][-1].rough.setp(2.9, vary=False)
        lipids[i].rough_t_a.constraint = structures[i][-1].rough
        lipids[i].rough_h_t.constraint = structures[i][-1].rough

    lipids, structures = two_layer.set_constraints(lipids, structures)

    models = toolbox.get_models(num_contrasts, structures, datasets)

    global_objective = toolbox.get_objectives(num_contrasts, models, datasets)

    if logl:
        to_return = global_objective.logl()
    else:
        to_return = global_objective
    return to_return
