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
        'paper/paper.pdf',
        'paper/si.pdf'

rule clear:
    run:
        shell('rm paper/*.pdf')
        shell('rm paper/*.txt')
        shell('rm output/*')

rule make_si:
    input:
        'paper/si.tex',
        'paper/ev_table.txt'
    output:
        'paper/si.pdf'
    run:
        shell("cd paper && pdflatex si.tex && pdflatex si.tex && cd ../")

rule make_paper:
    input:
        'paper/evidence.pdf',
        'paper/refl.pdf',
        'paper/iterations.pdf',
        'paper/post_all.pdf',
        'paper/post_best.pdf',
        'paper/best_per.pdf',
        ['paper/{}_ev.txt'.format(i) for i in VARIABLES],
        'paper/ev_table.txt',
        'paper/best_ev.txt',
        'paper/next_best_ev.txt',
        'paper/best_label.txt',
        'paper/next_best_label.txt',
        'paper/diff_ev.txt',
        'paper/d_h_range.txt',
        'paper/V_h_range.txt',
        'paper/d_t_range.txt',
        'paper/V_t_range.txt',
        'paper/paper.tex',
        'paper/paper.bib',
    output:
        'paper/paper.pdf'
    run:
        shell("cd paper && pdflatex paper.tex && bibtex paper.aux && pdflatex paper.tex && pdflatex paper.tex && cd ../")

rule evidence_plot:
    input:
        ['output/model^{}.h5'.format(i) for i in VARIABLES],
        'scripts/plotting.py'
    output:
        'paper/evidence.pdf',
        'paper/refl.pdf',
        'paper/iterations.pdf',
        'paper/post_all.pdf',
        'paper/post_best.pdf',
        'paper/best_per.pdf',
        ['paper/{}_ev.txt'.format(i) for i in VARIABLES],
        'paper/ev_table.txt',
        'paper/best_ev.txt',
        'paper/next_best_ev.txt',
        'paper/best_label.txt',
        'paper/next_best_label.txt',
        'paper/diff_ev.txt',
        'paper/d_h_range.txt',
        'paper/V_h_range.txt',
        'paper/d_t_range.txt',
        'paper/V_t_range.txt'
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