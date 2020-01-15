import json
import itertools

variables = [
    "th",
    "mvh",
    "tt", 
    "phit",
    "mvt",
    "rough", 
]

VARIABLES = []
for i in range(len(variables)):
    for j in itertools.combinations(variables, i + 1):
        VARIABLES.append('_'.join(j))

rule targets:
    input:
        ['output/model^{}.h5'.format(i) for i in VARIABLES]

rule analysis_all:
    input:
        'scripts/analysis.py',
        'scripts/two_layer.py',
        'scripts/toolbox.py',
        'scripts/model.py',
    output:
        ['output/model^{b}.h5'],
    run:
        shell("echo scripts/analysis.py {wildcards.b}")
        shell("python scripts/analysis.py {wildcards.b}")