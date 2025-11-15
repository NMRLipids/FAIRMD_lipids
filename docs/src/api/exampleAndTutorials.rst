
.. _exampleAndTutorials:

Examples and tutorials
======================

A template project for analyses utilizing the FAIRMD Lipids is available in `here <https://github.com/NMRLipids/databank-template/tree/main>`__ where examples and templates for analyses are available from the `scripts <https://github.com/NMRLipids/databank-template/tree/main/scripts>`__ folder.

#. `Plotting basic simulation properties <https://github.com/NMRLipids/databank-template/blob/main/scripts/plotSimulation.ipynb>`__

   Plots the basic properties of simulation selected based on its FAIRMD Lipids ID number. In addition to `GitHub <https://github.com/NMRLipids/databank-template/blob/main/scripts/plotSimulation.ipynb>`__, the code can be run at `Colab <https://colab.research.google.com/github/NMRLipids/databank-template/blob/main/scripts/plotSimulation.ipynb>`__ without additional setup, and app is available at `Streamlit <https://lolicato-nmrlipids-gui-app-qa2ffe.streamlit.app/>`__.

#. `Show ranking tables of simulations based in their quality against experimental data <https://github.com/NMRLipids/databank-template/blob/main/scripts/plotQuality.ipynb>`__

   Shows different kinds of rankings of simulations against experimental data. In addition to `GitHub <https://github.com/NMRLipids/databank-template/blob/main/scripts/plotQuality.ipynb>`__, the code can be run at `Colab <https://colab.research.google.com/github/NMRLipids/databank-template/blob/main/scripts/plotQuality.ipynb>`__ without additional setup, and app is available at `Streamlit <https://lolicato-nmrlipids-gui-app-qa2ffe.streamlit.app/>`__.

#. `Template for more advance API usage <https://github.com/NMRLipids/databank-template/blob/main/scripts/template.ipynb>`__

   Demonstrates the usage of API by three examples. 1) Selects a random simulation and prints the related databank content in human readable format. 2) Shows the readily analyzed properties for the selected random simulation (area per lipid, membrane thickness, relative equilibration times, X-ray scattering form factors, and C-H bond order parameters). 3) Selects a random simulation with the trajectory size below 100Mb and calculates P-N vector angle with respect to membrane normal for all lipids for which P and N atoms are available in headgroup. In addition to `GitHub <https://github.com/NMRLipids/databank-template/blob/main/scripts/template.ipynb>`__, the code can be run at `Colab <https://colab.research.google.com/github/NMRLipids/databank-template/blob/main/scripts/template.ipynb>`__ without additional setup.

#. `Showing statistics from the FAIRMD Lipids (OUTDATED) <https://github.com/NMRLipids/MuseumDatabank/blob/main/Scripts/AnalyzeDatabank/stats.ipynb>`__

   Plots distributions of simulation lengths, number of atoms and trajectory sizes, distribution of number of lipid components, available single component bilayers and binary lipid mixtures, and pie diagram of temperatures.
