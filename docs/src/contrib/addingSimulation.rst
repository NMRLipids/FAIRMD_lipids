.. _addSimulation:

Adding simulations
==================

Stepwise instructions to add simulation into the FAIRMD Lipids, run the basic analyses
and perform automatic quality evaluation are given here. The first three steps are a
minimum requirements to add a simulation. The first three steps can be performed using
graphical GitHub interface. To run the basic analyses and quality evaluation from steps
4 forward, you need to create a local fork of the `NMRlipids BilayerData git
<https://github.com/NMRLipids/BilayerData/>`_. 

#. Add trajectory and topology (tpr for Gromacs, pdb or corresponding to other programs)
   file into a `Zenodo <https://zenodo.org/>`_ repository.\ If you want to use other
   repository than Zenodo, please do not hesitate to open an `GitHub issue
   <https://github.com/NMRLipids/FAIRMD_lipids/issues>`_ on this.

#. Create an ``info.yaml`` file containing the essential information on your simulation
   by filling the `info template
   <https://github.com/NMRLipids/FAIRMD_lipids/blob/main/src/fairmd/lipids/SchemaValidateion/Schema/info_template.yaml>`_.
   For instructions, see :ref:`readmesimu` and `examples
   <https://github.com/NMRLipids/BilayerData/tree/main/info_files>`_. Mapping files are
   described in  :ref:`molecule_names` and are located in the :ref:`molecule_record`
   inside the folder of corresponding molecule.

#. You can store the created ``info.yaml`` file somewhere inside `./info_files/
   <https://github.com/NMRLipids/BilayerData/tree/main/info_files>`_ folder in the
   BilayerData git and make a pull request to the main branch. **You can stop here or
   continue to create ``README.yaml`` file in step 4.**

#. Before continuing, make sure that your BilayerData repository is switched to your own
   fork. To create the ``README.yaml`` file for the databank you should run
   :ref:`add_simulation_py` on your info-file:

   .. code-block:: bash

      fmdl_add_simulation -f {path to the info.yaml file that you created} -w
      {working-directory}

   After this is finished, you should see a new folder in ``./Simulations`` which contains
   the ``README.yaml`` file of your system. Commit the created ``README.yaml`` file into
   the git.

   .. code-block:: bash

      git add Simulations/**/README.yaml

#. To perform basic analyses for the added system(s) you should run
   :ref:`compute_databank_py` with the following keys:

   .. code-block:: bash

      fmdl_compute_databank --apl --nmrpca --ff --thickness \
         --OP --range *-0

   By default, your new-created system gets ID -1, -2, etc, so we run recomputing only
   on range from -inf to 0. After this, you should see the results in the same folder
   where ``README.yaml`` is located. The results files should be staged by running

   .. code-block:: bash

      git add Simulations/**/*.json

#. For the quality evaluation against experiments, the simulation needs to be first
   connected to the corresponding experimental data (if available) by running
   :ref:`match_experiments_py`. This will add the ``EXPERIMENT`` dictionary into the
   ``README.yaml`` files. This dictionary defines the location of related experimental
   data in `./experiments
   <https://github.com/NMRLipids/BilayerData/tree/main/experiments>`_ folder. Then the
   quality evaluation can be then done by running :ref:`quality_evaluation_py`:

   .. code-block:: bash

      fmdl_match_experiments fmdl_evaluate_quality

   The resulting qualities can be then added into the git by running

   .. code-block:: bash

      git add Simulations/**/README.yaml git add Simulations/**/*.json

   To create rankings of simulations based on their quality against experiments (see
   :ref:`make_ranking_py`) and to store the results in folder `Ranking
   <https://github.com/NMRLipids/BilayerData/tree/main/Ranking>`_, run

   .. code-block:: bash

      fmdl_make_ranking git add Ranking/*.json

#. Finally, commit the added data into your fork and make a pull request into the main
   branch.

.. toctree::
   :maxdepth: 1

   ../schemas/simulation_metadata.md
