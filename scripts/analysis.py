"""
Performs dynamic nested sampling for monolayer 

author: Andrew McCluskey (andrew.mccluskey@diamond.ac.uk)
"""

import sys
import numpy as np
import dynesty
import h5py
import model as nm

OUTPUT_DIR = "output/"

CON_LIST = ["d13acmw", "d13d2o", "hd2o", "d70acmw", "d70d2o", "d83acmw", "d83d2o"]
VARIABLES = sys.argv[1]
VAR_LIST = VARIABLES.split("_")

np.random.seed(1)


def ptform(uniform):
    """
    Prior informations for dynesty.

    Parameters
    ----------
    uniform : nd_array
        An array with length of number of parameters being varied.

    Returns
    -------
    prior : float
        Prior probability in a given range.
    """
    priors = []
    for i, v in enumerate(VAR_LIST):
        if "th" == v:
            priors.append(8 + (16 - 8) * uniform[i])
        if "mvh" == v:  
            priors.append(300 + (380 - 300) * uniform[i])
        if "tt" == v:  
            priors.append(10 + (26 - 10) * uniform[i])
        if "phit" == v:  
            priors.append(0.5 + (1 - 0.5) * uniform[i])
        if "mvt" == v:  
            priors.append(800 + (1000 - 800) * uniform[i])
        if "rough" == v:  
            priors.append(2.9 + (5 - 2.9) * uniform[i])
    return priors


def main():
    """
    Dynamic Nested Sampling of the Campbell model 3
    """
    file_name = "/model^{}".format(VARIABLES)

    ndim = len(VAR_LIST)

    variables = {}
    for i, v in enumerate(VAR_LIST):
        variables[v] = i

    sampler = dynesty.NestedSampler(
        nm.refl, ptform, ndim, logl_args=[CON_LIST, variables]
    )
    sampler.run_nested(dlogz=0.5, print_progress=True)
    res = sampler.results
    print(res)

    h5_file = h5py.File(OUTPUT_DIR + file_name + ".h5", "w")
    h5_file.create_dataset('samples', data=res['samples'])
    h5_file.create_dataset('logl', data=res['logl'])
    h5_file.create_dataset('niter', data=res['niter'])
    h5_file.create_dataset('logwt', data=res['logwt'])
    h5_file.create_dataset('logz', data=res['logz'])    
    h5_file.create_dataset('logzerr', data=res['logzerr'])    
    h5_file.create_dataset('logvol', data=res['logvol'])
    h5_file.close()


if __name__ == "__main__":
    main()
