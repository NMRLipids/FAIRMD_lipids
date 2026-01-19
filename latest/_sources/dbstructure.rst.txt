.. include:: global.rst

.. _dbstructure:

Data Organization
=================

The Databank includes three types of objects: Simulation, Experiment, and Molecule. Each is described by it's own metadata file. Metadata is stored in YAML format. Simulation object stores precomputed properties in JSON format. Experiment object stores experimental data in JSON format.

Databank objects are stored in the separated repository; by default, it's a `BilayerData <https://github.com/NMRlipids/BilayerData>`_ but any of its forks or even completely new repo could be used instead.
It is connected into the FAIRMD Lipids repository as a submodule into :file:`./Data` folder.

.. code-block:: console

   Data/
   ├── Simulations/
   │   └── rw5/
   │       └── g8v/
   │           └── au3/
   │               └── b8c/
   │                   ├── data-piece.json
   │                   └── README.yaml
   ├── experiments/
   │   ├── FormFactors/
   │   │   └── 10.1001/
   │   │       └── j.some.thing.3242.11/
   │   │           └── 1/
   │   │               ├── form-factor-data.json
   │   │               └── README.yaml
   │   └── OrderParameter/
   │   │   └── 10.1001/
   │   │       └── j.some.thing.3242.11/
   │   │           └── 1/
   │   │               ├── lipid1-data.json
   │   │               ├── lipid2-data.json
   │   │               └── README.yaml
   ├── info_files/
   │   └── some-folder/
   │       ├── info1.yaml
   │       └── info2.yaml
   ├── lipid_json_buildH/
   │   ├── ua-dictionary-1.json
   │   └── ua-dictionary-2.json
   ├── Ranking/
   │   ├── lipid1-headgroup-ranking.json
   │   └── lipid2-tail1-ranking.json
   └── Molecules/
       ├── membrane/
       │   └── lipid1/
       │       ├── lipid1-forcefield2-mapping.yaml
       │       └── metadata.yaml
       └── solution/
           └── ion1/
               ├── ion1-forcefield3-mapping.yaml
               └── metadata.yaml

.. _dbstructure_simu:

Simulation record
-----------------

Simulation metadata README.yaml
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
A file containing all the relevant information of each simulation entry. These files are
stored at :py:data:`fairmd.lipids.FMDL_SIMU_PATH`.

For more information see :doc:`schemas/simulation_metadata`.

Simulation computed properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. :file:`apl.json` contains area per lipid (Å\ :sup:`2`) as a function of time (ps) for
   a simulation. |br| See :py:func:`fairmd.lipids.analyze.computeAPL` for computing details.

2. :file:`XXXXOrderParameters.json` contains C-H bond order parameters and their
   uncertainties for each lipid XXXX.
   Key in json format is the merged universal C and H atom names. Three values are given:
   order parameter value for the C-H bond, standard deviation, and the standard error of
   the mean. |br| See :py:func:`fairmd.lipids.analyze.computeOP` for computing details.

3. :file:`FormFactor.json` contains simulated X-ray scattering form factor.
   X-axis value unit are Å\ :sup:`-1` and y-axis value units are e/nm\ :sup:`2`. |br| See
   :py:func:`fairmd.lipids.analyze.computeFF` for computing details.

4. :file:`XXXXDensity.json` contains electron densities calculated from trajectory.
   We store separately lipid, water, and total densities. X-axis units are nm and y-axis
   units e/nm\ :sup:`3`. |br| See :py:func:`fairmd.lipids.analyze.computeDensity` for
   computing details.

5. :file:`thickness.json` contains thickness (nm) calculated from the trajectory. |br|
   See :py:func:`fairmd.lipids.analyze.computeThickness` for computing details.

6. :file:`eq_times.json` contains a dictionary with relative equilibration time of
   molecule principal components. It demonstrates how long time is required in this
   particular simulation to sample the ensemble for each lipid molecule. |br|
   See :py:func:`fairmd.lipids.analyze.computeNMRPCA` for computing details.


.. _simulation-quality-files:
Simulation quality files
~~~~~~~~~~~~~~~~~~~~~~~~

1. :file:`XXXX_OrderParameters_quality.json` contains a dictionary with the quality of
   each C-H bond against experiments, if available, for molecule XXXX. First key is the
   DOI for the source of experimental data. Second key gives C-H pair univeral atom names.
   Third key gives values for simulation order parameter, its stardard deviation, standard
   error of the mean, experimental order parameter, its error, and finally the quality of
   the order parameter. Quality is the probability for the agreement between simulated and
   experimental values taking into account the error bars.

   For more details, see the `FAIRMD Lipids manuscript
   <https://doi.org/10.1038/s41467-024-45189-z>`_.

2. :file:`XXXX_FragmentQuality.json` contains fragment qualities determined separately
   for each lipid XXXX in the simulation with experimental data available using Eq. (4)
   in the `FAIRMD Lipids manuscript <https://doi.org/10.1038/s41467-024-45189-z>`_.

3. :file:`SYSTEM_quality.json` contains total qualities averaged over different lipids
   for different membrane parts calculated from Eq. (5) in the `FAIRMD Lipids manuscript
   <https://doi.org/10.1038/s41467-024-45189-z>`_.

4. :file:`FormFactorQuality.json` contains quality of form factor against experimental
   data determined as described in the `FAIRMD Lipids manuscript
   <https://doi.org/10.1038/s41467-024-45189-z>`_. Second term is the scaling
   coefficient for experimental intensities (Eq. (6) in the `FAIRMD Lipids manuscript
   <https://doi.org/10.1038/s41467-024-45189-z>`_).

.. toctree::
   :maxdepth: 1

   ../quality.rst


.. _dbstructure_exp:

Experiment record
-----------------

The data is located in subfolders named after the DOI of the original source publication
of the data. Because a publication can contain more than one experimental data set, each
data set is stored in a subfolder with integer name, for example,
`experiments/OrderParameters/10.1039/c2cp42738a/2
<https://github.com/NMRLipids/BilayerData/tree/main/experiments/OrderParameters/10.1039/c2cp42738a/2>`_.
Each of such folders then contain the experimental data and ``README.yaml`` files
describing the experimental conditions.

Solid-state NMR experiment metadata
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Order parameter (OP) values are measured by ssNMR, and the experiment folder contains
metadata file and measured OP values for each measured lipid.

1. :file:`README.yaml` contains metadata of the experiment, such as DOI, lamellar phase
   preparation protocol, temperature, B0 field, ions concentration and other relevant
   information (:ref:`explore all ssNMR metadata <readmeexp>`).

2. :file:`POPC_OrderParameters.json` contains C-H bond order parameters and their
   experimental uncertainties defined for universal atom names.

Form factor X-ray experiment metadata
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TODO: write the block

1. :file:`README.yaml` contains metadata of the experiment, such as DOI, properties of
   the source, type of liposome suspension, and other relevant information.

2. :file:`FormFactor.json` contains X-ray scattering form factor data.

.. _molecule_record:

Molecule record
---------------

Molecule metadata
~~~~~~~~~~~~~~~~~
Molecular metadata is stored in :file:`metadata.yaml` file, and its structure
is described in :ref:`molmetadata_schema`.

.. toctree::
   :maxdepth: 1
   :hidden:

   ../schemas/molecule_metadata.rst


Mapping files
~~~~~~~~~~~~~
Describe connection between universal atom names and simulation specific atom names.
These enable automatic analyses over simulations with different naming conventions. For
more information, see :ref:`molecule_names`

Other files
-----------

.. _ranking_files:

Ranking files
~~~~~~~~~~~~~~~~~
Contain tables where major conditions and quality scores are stored for each simulation in CSV format. These files are stored at ``Ranking`` subfolder of :py:data:`fairmd.lipids.FMDL_DATA_PATH` folder. We store:

1. Form-factor ranking files named as ``FF_ranking.csv``.
2. Order paparmeter ranking files named as ``OP_ranking.csv``.
3. Individual lipid ranking files named as ``<lipid>_ranking.csv``.

The difference between 2. and 3. is that 2. contains average quality scores over all lipids in the simulation, while 3. contains quality scores for only selected lipid.

.. _uadic_files:

Lipid united-atom dictionary
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Contains information about guessing implicit hydrogens for united-atom force fields. See :ref:`NMRlipids/buildh <https://github.com/NMRlipids/buildh>`_ for more details.

TODO: write the block

.. _info_files:

info.yaml
~~~~~~~~~
Contains information given by contributor when adding data into the FAIRMD Lipids, i.e.
all non-recomputable fields of simulation ``README.yaml``. This file is given as an
input to :ref:`add_simulation_py` to create ``README.yaml`` files:
:code:`fmdl_add_simulation -f info.yaml`. These files are not required but currently are
stored historically in ``info_files`` subfolder of
:py:data:`fairmd.lipids.FMDL_DATA_PATH`.

For more information see :ref:`readmesimu` and :ref:`addSimulation`.

