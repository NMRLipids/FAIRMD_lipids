# FAIRMD Lipids

[![tests status](https://img.shields.io/github/checks-status/NMRLipids/FAIRMD_lipids/main)](https://github.com/NMRLipids/FAIRMD_lipids/actions?query=branch%3Amain)
[![documentation-stable](https://img.shields.io/badge/📚_documentation-stable-sucess)](https://nmrlipids.github.io/FAIRMD_lipids/)
[![documentation-latest](https://img.shields.io/badge/📒_documentation-latest-yellow)](https://nmrlipids.github.io/FAIRMD_lipids/latest/index.html)
[![coverage](https://codecov.io/gh/NMRLipids/FAIRMD_lipids/branch/main/graph/badge.svg)](https://codecov.io/gh/NMRLipids/FAIRMD_lipids)

This is the FAIRMD Lipids &mdash; an API package for the community-driven catalogue 
FAIRMD Lipid Databank (formerly NMRlipids Databank) containing atomistic MD
simulations of biologically relevant lipid membranes emerging from the [NMRlipids Open
Collaboration](http://nmrlipids.blogspot.com/2021/03/second-online-meeting-on-nmrlipids.html). 

# Installation

The code has been tested in Linux and MacOS environment with python 3.10 or newer.
Recent [Gromacs](https://manual.gromacs.org/current/install-guide/index.html) version
should be available in the system. All dependecies are listed in
[pyproject.toml](pyproject.toml).

We recommend installing python libraries into an environment, for example, using conda:

```bash
conda create --name fairmd-lipids python==3.10 -c conda-forge
conda activate fairmd-lipids
```

Install *fairmd-lipids* package from repo:

```bash
pip install git+https://github.com/NMRLipids/FAIRMD_lipids.git
```

or from pypi:

```bash
pip install fairmd-lipids
```

Note that the data is stored as a separated repository and should be loaded after
cloning. Default data storage is
[BilayerData](https://github.com/NMRLipids/BilayerData). You **MUST** specify
`FMDL_DATA_PATH` before start working. The easiest way to start is to use
`fmdl_initialize_data` script provided with the package:

```bash
fmdl_initialize_data stable
source databank_env.rc
```

Then you can work with the standalone scripts as well as use `fairmd-lipids` package in
your python code.

# Documentation

The FAIRMD Lipids project documentation is available in
[here](https://nmrlipids.github.io/FAIRMD_lipids/latest/). More information and example
applications are available from the [FAIRMD Lipids
manuscript](https://doi.org/10.1038/s41467-024-45189-z).

The `fairmd-lipids` python module provides programmatic access to all simulation data in
the FAIRMD Lipids. It allows to request data to construct various datasets on the base
of experimental and simulation data from the database that allow one to learn various
models about lipid bilayer properties. It also allows to design and effectively run
automated analysis across all the simulations in the database.

## How to use 

A [jupyter template
notebook](https://github.com/NMRLipids/databank-template/blob/main/scripts/template.ipynb)
can be used to get started with the analyses utilizing the FAIRMD Lipids.

Connection of [Universal molecule and atom naming
conventions](https://nmrlipids.github.io/FAIRMD_lipids/latest/schemas/moleculesAndMapping.html)
with simulation specific names delivered by mapping files can be used to perform
automatic analyses over large sets of simulations. The results for large analyses can be
stored using the same structure as used for `README.yaml` files as done, for example,
for [water
permeation](https://github.com/NMRLipids/FAIRMD_lipids/tree/main/Data/MD-PERMEATION)
and lipid
[flip-flop](https://github.com/NMRLipids/DataBankManuscript/tree/main/Data/Flipflops)
rates in the [repository related to the FAIRMD Lipids
manuscript](https://github.com/NMRLipids/DataBankManuscript).

# Web UI

[FAIRMD Lipids-webUI](https://databank.nmrlipids.fi/) provides an easy access to the
FAIRMD Lipids content. Simulations can be searched based on their molecular composition,
force field, temperature, membrane properties, and quality; the search results are
ranked based on the simulation quality as evaluated against experimental data when
available. Web-UI provides basic graphical reports for the computed properties as well
as graphical comparison between simulation and experimental data.

The Web-UI is being developed in the repository
[BilayerGUI_laravel](https://github.com/NMRlipids/BilayerGUI_laravel).

# Contribution

The project is open for contributions! 

Please consult [CONTRIBUTION.md](./CONTRIBUTION.md) for further information.

Do not hesitate contacting us via [databank@nmrlipids.fi](email:databank@nmrlipids.fi) if 
you have any questions or inquiries!
