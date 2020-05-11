# A general approach to maximise information density in neutron reflectometry analysis

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3697795.svg)](https://doi.org/10.5281/zenodo.3697795) [![arXiv](https://img.shields.io/badge/arXiv-1910.10581-orange.svg)](https://arxiv.org/abs/1910.10581) [![License](https://img.shields.io/github/license/arm61/model_select.svg?color=lightgrey)](https://github.com/arm61/model_select/blob/master/LICENSE)

[Andrew R. McCluskey](https://orcid.org/0000-0003-3381-5911)&ast;, [Thomas Arnold](https://orcid.org/0000-0001-7196-7831), Joshaniel F. K. Cooper, and [Tim Snow](https://orcid.org/0000-0001-7146-6885)

&ast;[andrew.mccluskey@diamond.ac.uk](mailto:andrew.mccluskey@diamond.ac.uk)/[a.r.mccluskey@bath.ac.uk](mailto:a.r.mccluskey@bath.ac.uk)

*Outlining and applying a Bayesian model selection framework for neutron reflectometry analysis.*

This is the electronic supplementary information (ESI) associated with the publication "Using Bayesian model selection to advise neutron reflectometry analysis from Langmuir-Blodgett monolayers".
This ESI provides full details of the analyses performed in the work and access to an automated analysis workflow, through this we aim to provide better access to analysis reproducibility.
For more information about reproducible workflows in Python, check out [Tania Allard's talk from Pycon2018](http://bitsandchips.me/Talks/PyCon.html#/title).

## [Data](./data)

The reduced neutron reflectometry data can be found in this repository, in the [data](./data) directory.

Note, that the data was originally collected by Hollinshead and co-workers (DOI: [10.1021/la8028319](https://doi.org/10.1021/la8028319)).

## Analysis

This ESI aims to provide a fully reproducible workflow to the data analysis presented within the paper.

Requirements:

- anaconda or miniconda python
- [REVTeX](https://journals.aps.org/revtex)

The supplied Snakefile, will reproduce all of the analysis, plot the figures, and build a preprint version of the paper (`paper/paper.pdf`) when run.
Be aware that the analyses within this work are non-trivial and take many hours to run so **use caution** before re-running.

To re-run all analysis, simply run the following commands:

```
conda env create --prefix ./model_select --file enviroment.yml

source activate ./model_select

snakemake clear # this will remove all of the output from previous runs

snakemake
```

The [`Snakefile`](./Snakefile) outlines the process of this reproducible analysis and [all code](./scripts) used is available. 

## [Figures](./paper)

PDF versions of the figures, can be found in the [`paper`](./paper) directory.--->

## Acknowledgements

A.R.M. would like to acknowledge David J. Barlow and M. Jayne Lawrence for kindly sharing the neutron reflectometry data used, and Simon Titmuss for suggesting the use of the `dynesty` package.
This work is supported by the Ada Lovelace Centre – a joint initiative between the Science and Technology Facilities Council (as part of UK Research and Innovation), Diamond Light Source, and the UK Atomic Energy Authority.
