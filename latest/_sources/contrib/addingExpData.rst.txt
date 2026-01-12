.. _addingExpData:

Adding experimental data
========================

Experimental data is stored in :ref:`BilayerData/experiments <dbstructure_exp>` folder. C-H bond order parameters from NMR are in ``OrderParameters`` subfolder and X-ray scattering form factors in ``FormFactors`` subfolder.
The keys of these dictionaries are summarized in the :ref:`Experiment metadata description <readmeexp>`.

**Steps to add experimental data**

#. Create and fill the ``README.yaml`` file of your data.

#. Copy this README file data into a appropriate directory named as described above.

#. If you have order parameter data, create a file named ``{lipidname}_OrderParameters.json``
   where ``{lipidname}`` is the universal name of the lipid from which the data is measured
   from. The first two columns of this file should define the atom pair with universal
   atom names, third column has the experimental order parameter value, and fourth
   column has the experimental error. If the experimental error is not known, set it to 0.02.
   Store the created ``{lipidname}_OrderParameters.json`` file into the appropriate folder
   with the ``README.yaml`` file. Create ``json`` version from ``dat`` file by running

   .. code-block:: bash

        python data_to_json.py path-to-dat.dat

   in folder `experiments/OrderParameters <https://github.com/NMRLipids/BilayerData/tree/main/experiments/OrderParameters>`_. You can see previously added experiments for examples.

#. If you have X-ray scattering form factor data, store the form factor into appropriate
   folders in ASCII format where first column in x-axis values (Å\ :sup:`-1`), second
   column is y-axis value, and third is the error. Then create ``json`` file by runnig

   .. code-block:: bash

      python data_to_json.py ascii-file.txt

   in folder [experiments/FormFactors](https://github.com/NMRLipids/BilayerData/tree/main/experiments/FormFactors).
   Please see previously added experiments for examples.

#. Adding experiments can lead to recalculation of quality of some simulations and rankings.
   So, after the addition, please run

   .. code-block:: bash

      fmdl_match_experiments
      fmdl_evaluate_quality
      fmdl_make_ranking

#. Submit the files to your branch and make a pull-request.

.. toctree::
   :maxdepth: 1

   ../schemas/experiment_metadata.md

