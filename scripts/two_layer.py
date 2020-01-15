"""
OneLayer class

author: Andrew McCluskey (andrew.mccluskey@diamond.ac.uk)
"""
import numpy as np
from refnx.analysis import possibly_create_parameter, Parameters
from refnx.reflect import Component
import toolbox

# pylint: disable=R0902


def get_b_lipid(contrasts):
    """
    Get the scattering length for the lipid for given contrasts

    Parmaeters
    ----------
    contrasts : list
        List of contrasts

    Returns
    -------
    b_lipid : list of complex
        List of scattering lengths for the given component in each contrast
    """
    b_lipid = []
    for cont in contrasts:
        if cont[:3] == "d13":
            b_lipid.append(
                [
                    toolbox.get_scattering_length(
                        {"C": 10, "D": 13, "H": 5, "O": 8, "N": 1, "P": 1}
                    ),
                    toolbox.get_scattering_length(
                        {"C": 17 * 2, "H": 17 * 4 + 2}
                    ),
                ]
            )
        elif cont[:3] == "d70":
            b_lipid.append(
                [
                    toolbox.get_scattering_length(
                        {"C": 10, "D": 5, "H": 13, "O": 8, "N": 1, "P": 1}
                    ),
                    toolbox.get_scattering_length(
                        {"C": 17 * 2, "D": 17 * 4 + 2}
                    ),
                ]
            )
        elif cont[:3] == "d83":
            b_lipid.append(
                [
                    toolbox.get_scattering_length(
                        {"C": 10, "D": 18, "O": 8, "N": 1, "P": 1}
                    ),
                    toolbox.get_scattering_length(
                        {"C": 17 * 2, "D": 17 * 4 + 2}
                    ),
                ]
            )
        elif cont[:1] == "h":
            b_lipid.append(
                [
                    toolbox.get_scattering_length(
                        {"C": 10, "H": 18, "O": 8, "N": 1, "P": 1}
                    ),
                    toolbox.get_scattering_length(
                        {"C": 17 * 2, "H": 17 * 4 + 2}
                    ),
                ]
            )
    return b_lipid


class TwoLayer(Component):
    """
    The class to describe a two-layer model for the lipid monolayer.

    Parameters
    ----------
    bs: float, array_like
        The scattering lengths for the head and tail components
    name: string
        A name for the monolayer
    """

    def __init__(self, bs, name="two_layer"):
        super(TwoLayer, self).__init__()
        if isinstance(bs[0], complex):
            self.b_real_h = possibly_create_parameter(
                bs[0].real, "{} - b_real_head".format(name)
            )
            self.b_imag_h = possibly_create_parameter(
                bs[0].imag, "{} - b_imag_head".format(name)
            )
        else:
            self.b_real_h = possibly_create_parameter(
                bs[0], "{} - b_real_head".format(name)
            )
            self.b_imag_h = possibly_create_parameter(
                0, "{} - b_imag_head".format(name)
            )
        if isinstance(bs[1], complex):
            self.b_real_t = possibly_create_parameter(
                bs[1].real, "{} - b_real_tail".format(name)
            )
            self.b_imag_t = possibly_create_parameter(
                bs[1].imag, "{} - b_imag_tail".format(name)
            )
        else:
            self.b_real_t = possibly_create_parameter(
                bs[1], "{} - b_real_tail".format(name)
            )
            self.b_imag_t = possibly_create_parameter(
                0, "{} - b_imag_tail".format(name)
            )

        self.mol_vol_h = possibly_create_parameter(
            100, "{} - molecular_volume_head".format(name)
        )
        self.mol_vol_t = possibly_create_parameter(
            100, "{} - molecular_volume_tail".format(name)
        )

        self.thick_h = possibly_create_parameter(
            100, "{} - thickness_head".format(name)
        )
        self.thick_t = possibly_create_parameter(
            100, "{} - thickness_tail".format(name)
        )

        self.phi_h = possibly_create_parameter(
            0.5, "{} - solvation_head".format(name)
        )
        self.phi_t = possibly_create_parameter(
            0, "{} - solvation_tail".format(name)
        )

        self.rough_h_t = possibly_create_parameter(
            3.3, "{} - roughness_head_tail".format(name)
        )
        self.rough_t_a = possibly_create_parameter(
            3.3, "{} - roughness_tail_air".format(name)
        )

        self.name = name

    def slabs(self, structure=None):
        """
        Returns
        -------
        slab_model = array of np.ndarray
            Slab representation of monolayer made up of two layers
        """
        layers = np.zeros((2, 5))

        layers[0, 0] = self.thick_t
        layers[0, 1] = self.b_real_t * 1e6 / self.mol_vol_t
        layers[0, 2] = self.b_imag_t * 1e6 / self.mol_vol_t
        layers[0, 3] = self.rough_t_a
        layers[0, 4] = 1 - self.phi_t

        layers[1, 0] = self.thick_h
        layers[1, 1] = self.b_real_h * 1e6 / self.mol_vol_h
        layers[1, 2] = self.b_imag_h * 1e6 / self.mol_vol_h
        layers[1, 3] = self.rough_h_t
        layers[1, 4] = 1 - self.phi_h

        return layers

    @property
    def parameters(self):
        para = Parameters(name=self.name)
        para.extend(
            [
                self.b_real_h,
                self.b_imag_h,
                self.b_real_t,
                self.b_imag_t,
                self.thick_h,
                self.thick_t,
                self.mol_vol_h,
                self.mol_vol_t,
                self.rough_h_t,
                self.rough_t_a,
                self.phi_h,
                self.phi_t,
            ]
        )
        return para

    def logp(self):
        if self.phi_h >= 1 or self.phi_h < 0:
            return -np.inf
        if self.phi_t >= 1 or self.phi_t < 0:
            return -np.inf
        return 0


def set_constraints(lipids, structures):
    """
    Assign necessary constraints

    Parameters
    ----------
    lipids : OneLayer
        Lipids objects to constrained
    structures : refnx.reflect.Structure
        Structures to be constrained

    Returns
    -------
    lipids : OneLayer
        Constrained lipids objects
    structures : refnx.reflect.Structure
        Constrained structures
    """
    for i in range(1, len(lipids)):
        lipids[i].thick_h.constraint = lipids[0].thick_h
        lipids[i].thick_t.constraint = lipids[0].thick_t
        lipids[i].mol_vol_h.constraint = lipids[0].mol_vol_h
        lipids[i].mol_vol_t.constraint = lipids[0].mol_vol_t
        lipids[i].rough_h_t.constraint = lipids[0].rough_h_t
        lipids[i].rough_t_a.constraint = lipids[0].rough_t_a
        lipids[i].phi_h.constraint = lipids[0].phi_h
        lipids[i].phi_t.constraint = lipids[0].phi_t
    for i in range(1, len(structures)):
        structures[i][-1].rough.constraint = structures[0][-1].rough
    return lipids, structures
