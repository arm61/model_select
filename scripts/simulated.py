"""
Using nested sampling on simulated layer data.

Copyright (c) Andrew R. McCluskey

Distributed under the terms of the MIT License

@author: Andrew R. McCluskey
"""

# pylint: disable=R0903

import sys
import h5py
import numpy as np
import dynesty
from refnx.reflect import ReflectModel, SLD
from refnx.dataset import Data1D
from refnx.analysis import Transform, Objective
from multiprocessing import Pool, cpu_count

np.random.seed(1)


class SimulatedData:
    """
    Data simulation for two layers.
    """
    def __init__(self):
        """
        Initialisation function.
        """
        self.q_vectors = np.linspace(0.05, 0.5, 100)
        air = SLD(0, 'air')
        water = SLD(6.35, 'D2O')
        layer_1 = SLD(2., 'layer_1')
        layer_2 = SLD(4., 'layer_2')

        self.model = air(0, 0) | layer_1(25, 0) | layer_2(25, 0) | water(0, 0)

    def get_refl(self, uncertainty_mod=1e-13):
        """
        Return the reflectometry.

        Args:
            uncertainty_mod (float): Modulation of uncertainty.

        Returns:
            (array_like) Simulated reflected intensity.
        """
        refl_model = ReflectModel(self.model)
        refl_model.scale.setp(1, vary=False)
        refl_model.bkg.setp(5e-7, vary=False)
        return (
            self.q_vectors,
            refl_model(self.q_vectors),
            uncertainty_mod/refl_model(self.q_vectors),
        )


def structure_one(values):
    """
    Structure with one layer.

    Args:
        values (list): list of parameter values.

    Returns:
        (refnx.reflect.Structure) The structure with one layer.
    """
    air = SLD(0, 'air')
    water = SLD(6.35, 'D2O')
    layer_1 = SLD(values[0], 'layer_1')

    layer_1_layer = layer_1(values[1], 0)

    layer_1_layer.thick.setp(vary=False)
    layer_1_layer.sld.real.setp(vary=False)
    layer_1_layer.rough.setp(vary=False)

    structure = air | layer_1_layer | water
    return structure


def structure_two(values):
    """
    Structure with two layers.

    Args:
        values (list): list of parameter values.

    Returns:
        (refnx.reflect.Structure) The structure with two layers.
    """
    air = SLD(0, 'air')
    water = SLD(6.35, 'D2O')
    layer_1 = SLD(values[0], 'layer_1')
    layer_2 = SLD(values[1], 'layer_2')

    layer_1_layer = layer_1(values[2], 0)
    layer_2_layer = layer_2(values[3], 0)

    layer_1_layer.thick.setp(vary=False)
    layer_1_layer.sld.real.setp(vary=False)
    layer_1_layer.rough.setp(vary=False)

    layer_2_layer.thick.setp(vary=False)
    layer_2_layer.sld.real.setp(vary=False)
    layer_2_layer.rough.setp(vary=False)

    structure = air | layer_1_layer | layer_2_layer | water
    return structure


def structure_three(values):
    """
    Structure with three layers.

    Args:
        values (list): list of parameter values.

    Returns:
        (refnx.reflect.Structure) The structure with three layers.
    """
    air = SLD(0, 'air')
    water = SLD(6.35, 'D2O')
    layer_1 = SLD(values[0], 'layer_1')
    layer_2 = SLD(values[1], 'layer_2')
    layer_3 = SLD(values[2], 'layer_3')

    layer_1_layer = layer_1(values[3], 0)
    layer_2_layer = layer_2(values[4], 0)
    layer_3_layer = layer_3(values[5], 0)

    layer_1_layer.thick.setp(vary=False)
    layer_1_layer.sld.real.setp(vary=False)
    layer_1_layer.rough.setp(vary=False)

    layer_2_layer.thick.setp(vary=False)
    layer_2_layer.sld.real.setp(vary=False)
    layer_2_layer.rough.setp(vary=False)

    layer_3_layer.thick.setp(vary=False)
    layer_3_layer.sld.real.setp(vary=False)
    layer_3_layer.rough.setp(vary=False)

    structure = air | layer_1_layer | layer_2_layer | layer_3_layer | water
    return structure


def structure_four(values):
    """
    Structure with four layers.

    Args:
        values (list): list of parameter values.

    Returns:
        (refnx.reflect.Structure) The structure with four layers.
    """
    air = SLD(0, 'air')
    water = SLD(6.35, 'D2O')
    layer_1 = SLD(values[0], 'layer_1')
    layer_2 = SLD(values[1], 'layer_2')
    layer_3 = SLD(values[2], 'layer_3')
    layer_4 = SLD(values[3], 'layer_4')

    layer_1_layer = layer_1(values[4], 0)
    layer_2_layer = layer_2(values[5], 0)
    layer_3_layer = layer_3(values[6], 0)
    layer_4_layer = layer_4(values[7], 0)

    layer_1_layer.thick.setp(vary=False)
    layer_1_layer.sld.real.setp(vary=False)
    layer_1_layer.rough.setp(vary=False)

    layer_2_layer.thick.setp(vary=False)
    layer_2_layer.sld.real.setp(vary=False)
    layer_2_layer.rough.setp(vary=False)

    layer_3_layer.thick.setp(vary=False)
    layer_3_layer.sld.real.setp(vary=False)
    layer_3_layer.rough.setp(vary=False)

    layer_4_layer.thick.setp(vary=False)
    layer_4_layer.sld.real.setp(vary=False)
    layer_4_layer.rough.setp(vary=False)

    structure = (
        air | layer_1_layer | layer_2_layer |
        layer_3_layer | layer_4_layer | water
    )
    return structure


def analysis(values, n_layers, logl=True):
    """
    Evalutate the reflectometry from two-layer phospholipid monolayer model.

    Args:
        values (array_like): An array of values for the parameters that are
            being varied.
        logl (bool, optional): Return just the logl value (if `True`) or the
            whole refnx.analysis.GlobalObjective object (if `False`).

    Returns:
        (float or refnx.analysis.GlobalObjective) Depending on the value of
            `logl`, either the logl of the input or the whole
            GlobalObjective is returned.
    """
    if n_layers == 1:
        structure = structure_one
    elif n_layers == 2:
        structure = structure_two
    elif n_layers == 3:
        structure = structure_three
    elif n_layers == 4:
        structure = structure_four
    else:
        raise NotImplementedError("Max number of layers is 4")

    model = ReflectModel(structure(values), dq=5.0)
    model.scale.setp(1, vary=False)
    model.bkg.setp(5e-7, vary=False)

    sim_data = SimulatedData()

    data = Data1D(data=(sim_data.get_refl()))
    objective = Objective(model, data, transform=Transform('logY'))
    if logl:
        return objective.logl()
    return objective


def ptform_one(uniform):
    """
    Prior informations for dynesty with two layers.

    Args:
        uniform (np.nd_array): An array with length of number of parameters
            being varied.

    Returns:
        prior (list): Prior probability for each of the varying parameters.
    """
    priors = []
    priors.append(6.35 * uniform[0])
    priors.append(50 * uniform[1])
    return priors


def ptform_two(uniform):
    """
    Prior informations for dynesty with two layers.

    Args:
        uniform (np.nd_array): An array with length of number of parameters
            being varied.

    Returns:
        prior (list): Prior probability for each of the varying parameters.
    """
    priors = []
    priors.append(6.35 * uniform[0])
    priors.append(6.35 * uniform[1])
    priors.append(50 * uniform[2])
    priors.append(50 * uniform[3])
    return priors


def ptform_three(uniform):
    """
    Prior informations for dynesty with three layers.

    Args:
        uniform (np.nd_array): An array with length of number of parameters
            being varied.

    Returns:
        prior (list): Prior probability for each of the varying parameters.
    """
    priors = []
    priors.append(6.35 * uniform[0])
    priors.append(6.35 * uniform[1])
    priors.append(6.35 * uniform[2])
    priors.append(50 * uniform[3])
    priors.append(50 * uniform[4])
    priors.append(50 * uniform[5])
    return priors


def ptform_four(uniform):
    """
    Prior informations for dynesty with four layers.

    Args:
        uniform (np.nd_array): An array with length of number of parameters
            being varied.

    Returns:
        prior (list): Prior probability for each of the varying parameters.
    """
    priors = []
    priors.append(1 + (6.35 - 1) * uniform[0])
    priors.append(1 + (6.35 - 1) * uniform[1])
    priors.append(1 + (6.35 - 1) * uniform[2])
    priors.append(1 + (6.35 - 1) * uniform[3])
    priors.append(10 + (50 - 10) * uniform[4])
    priors.append(10 + (50 - 10) * uniform[5])
    priors.append(10 + (50 - 10) * uniform[6])
    priors.append(10 + (50 - 10) * uniform[7])
    return priors


def main(n_layers):
    """
    Nested Sampling of simple simulated data.

    Args:
        n_layers (int): Number of layers in analysis.
    """
    file_name = "output/simulated^{}".format(n_layers)

    ndim = n_layers * 2

    if n_layers == 1:
        ptform = ptform_one
    elif n_layers == 2:
        ptform = ptform_two
    elif n_layers == 3:
        ptform = ptform_three
    elif n_layers == 4:
        ptform = ptform_four
    else:
        raise NotImplementedError("Max number of layers is 4")
    p = Pool(cpu_count())
    print("CPUS", cpu_count())
    sampler = dynesty.DynamicNestedSampler(
        analysis, ptform, ndim, logl_args=[n_layers], pool=p, queue_size=cpu_count(),
    )
    sampler.run_nested(print_progress=True, wt_kwargs={'pfrac': 0.0}, stop_kwargs={'pfrac': 0.0})
    res = sampler.results
    print(res)

    h5_f = h5py.File(file_name + ".h5", "w")
    h5_f.create_dataset('samples', data=res['samples'])
    h5_f.create_dataset('logl', data=res['logl'])
    h5_f.create_dataset('niter', data=res['niter'])
    h5_f.create_dataset('logwt', data=res['logwt'])
    h5_f.create_dataset('logz', data=res['logz'])
    h5_f.create_dataset('logzerr', data=res['logzerr'])
    h5_f.create_dataset('logvol', data=res['logvol'])
    h5_f.close()


if __name__ == "__main__":
    main(int(sys.argv[1]))
