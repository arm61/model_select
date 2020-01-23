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
        'figures/evidence.pdf',
        'figures/refl.pdf',
        'figures/iterations.pdf',
        'figures/post_all.pdf',
        'figures/post_best.pdf',
        'figures/best_per.pdf',
        ['results/{}_ev.txt'.format(i) for i in VARIABLES],
        'results/ev_table.txt',
        'results/best_ev.txt',
        'results/next_best_ev.txt',
        'results/best_label.txt',
        'results/next_best_label.txt',
        'results/diff_ev.txt',
        'results/d_h_range.txt',
        'results/V_h_range.txt',
        'results/d_t_range.txt',
        'results/V_t_range.txt',
        ['output/simulated^{}.h5'.format(i) for i in range(1, 5)]

rule evidence_plot:
    input:
        ['output/model^{}.h5'.format(i) for i in VARIABLES],
        'scripts/plotting.py'
    output:
        'figures/evidence.pdf',
        'figures/refl.pdf',
        'figures/iterations.pdf',
        'figures/post_all.pdf',
        'figures/post_best.pdf',
        'figures/best_per.pdf',
        ['results/{}_ev.txt'.format(i) for i in VARIABLES],
        'results/ev_table.txt',
        'results/best_ev.txt',
        'results/next_best_ev.txt',
        'results/best_label.txt',
        'results/next_best_label.txt',
        'results/diff_ev.txt',
        'results/d_h_range.txt',
        'results/V_h_range.txt',
        'results/d_t_range.txt',
        'results/V_t_range.txt'
    run:
        shell("echo scripts/plotting.py")
        shell("python scripts/plotting.py")

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

rule analysis_simulated:
    input:
        'scripts/simulated.py',
    output:
        ['output/simulated^{b}.h5'],
    run:
        shell("echo scripts/analysis.py {wildcards.b}")
        shell("python scripts/analysis.py {wildcards.b}")