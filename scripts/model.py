"""
Model 3 of the Campbell models

author: Andrew McCluskey (andrew.mccluskey@diamond.ac.uk)
"""

from refnx.dataset import ReflectDataset
import toolbox
import two_layer


DATA_DIR = "data"


def refl(values, contrasts, variables, logl=True):
    """
    Evalutate the reflectometry from model 3

    Parameters
    ----------
    t_t : float
        Lipid tail thickness
    contrast : list
        Contrasts to be analysed
    logl: bool
        Should the logl likelihood be return, else the global objective is
        returned

    Returns
    -------
    if logl:
        logl : float
            Log-likelihood for the value of phi
    else:
        global_objective : refnx.analysis.GlobalObjective
            GlobalObjective for system
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
            structures[i][-1].rough.setp(values[variables["rough"]], vary=False)
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
