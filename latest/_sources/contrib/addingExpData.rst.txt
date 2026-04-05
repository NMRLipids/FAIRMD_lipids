.. _addingExpData:

Adding experimental data
========================

Experimental data is stored in :ref:`BilayerData/experiments <dbstructure_exp>` folder. C-H bond
order parameters from NMR are in ``OrderParameters`` subfolder and X-ray scattering form factors in
``FormFactors`` subfolder. The keys of these dictionaries are summarized in the :ref:`Experiment
metadata description <readmeexp>`.

**Additional material for the chapter:**

.. toctree::
   :maxdepth: 1

   ../schemas/experiment_metadata.md


Steps to add experimental data
------------------------------

#. Fork the BilayerData repository in the GitHub interface and add+checkout the branch with an explicit
   naming, e.g., ``add-exp-Smidth-2023-dpps``. Work only in this branch. You can add several experiments into
   one branch, it's also fine.

#. Create and fill the ``README.yaml`` file of your data. This is not as simple as it seems,
   please read the :ref:`experiment metadata guidelines <addexp_guidelines>` if you are not feeling familiar.

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

   in folder `experiments/OrderParameters
   <https://github.com/NMRLipids/BilayerData/tree/main/experiments/OrderParameters>`_. You can see
   previously added experiments for examples.

#. If you have X-ray scattering form factor data, store the form factor into appropriate
   folders in ASCII format where first column in x-axis values (Å\ :sup:`-1`), second
   column is y-axis value, and third is the error. Then create ``json`` file by runnig

   .. code-block:: bash

      python data_to_json.py ascii-file.txt

   in folder `experiments/FormFactors
   <https://github.com/NMRLipids/BilayerData/tree/main/experiments/FormFactors>`_. Please see
   previously added experiments for examples.

#. Adding experiments can lead to recalculation of quality of some simulations and rankings.
   So, after the addition, you may want to run

   .. code-block:: bash

      fmdl_match_experiments
      fmdl_evaluate_quality
      fmdl_make_ranking

   This is not obligatory. We will do it anyway at the CI/CD stage; however, it's useful to check that
   your experiments are properly paired and quality is recalculated.

#. Submit the files to your branch and make a pull-request.

.. _addexp_guidelines:

Guidelines to fill experiment metadata
--------------------------------------

Every field of the metadata file is :ref:`explained here <readmeexp>`. Our experience indicates that
reading the original paper and filling in the fields is not very straightforward, because
experimenters often use very different names to name the same things. Here, we will go through the
most strange points.

#. **Hydration**

   We use water mass % for hydration, i.e., 10 mg in 100 ml water is ``10/(10+100)*100`` (%).
   Experimentators could use whatever - v/w %, lipid:water ratio, or molar concentration. They should
   be converted to water mass %.

#. **Membrane composition**

   We use molar fractions for membrane composition. It should not be mixed with mass fraction. E.g. if
   the experiment has 1:1 (mass) TOCL and POPC, which means that it has 1:2 (mol) because TOCL is
   almost two times harder than POPC.

#. **Solution composition**

   We use mass fraction for solution composition, and we count each ion separately. If the membrane
   is neutral, the solution composition should be converted to mass % and written. For example, we
   have 150 mM NaCl, which means that you have ``0.150 * 35.5 = 5.33`` (g) in 1L and ``0.150 * 23 =
   3.5`` (g) in 1L. If it’s a diluted solution, then 1L is ~1000g, so we have 0.53% Cl- and 0.25%
   Na+.

   .. code-block:: yaml

      REAGENT_SOURCES:
         CLA: 0.53
         SOD: 0.25

   Life becomes more complex if you have an anionic lipid. Let’s imagine you have DOPS (sodium salt)
   added to 150 mM sodium chloride solution so that you have 50 mass % lipid. For a 1000 g lipid salt,
   you have 1000 g buffer. Na(DOPS) has MW=810 (g/mol), i.e., mass of Na in this salt is ``23/810*1000=28.4`` (g).
   So, in total, you have ``3.5 + 28.4`` (g) of Na+, i.e., the total amount of sodium is 3.2% and not 0.53%
   like in previous example. If the lipid solution is diluted, which is quite often in SAXS
   experiments, then the amount of cations from lipid salts is small. However, for solid-state NMR
   experiments, lipid concentration is always very high, as well as the abundance of counterions coming
   with them.

#. **Additional molecules**

   In the ideal world, all that we have in the sample vial, we should also have in the simulation box.
   However, this goal is not achievable. Some samples have additives, which are hard to simulate
   because of a lack of force fields or because they are too diluted to be added to a small simulation
   box. Molecules, which are hard to simulate but exist in the experiment, come into the category
   ``ADDITIONAL_MOLECULES``. Typical examples are buffers, antioxidants, chelators like EDTA, etc.

#. **Reagent sources**

   We care about reagent sources. Lipids can be bought synthetic or isolated from crude lipid extracts.
   They can be bought isolated, but are almost pure. Lipids can be bought partly deuterated and can be
   synthesized locally in the paper from synthetic lipids. Note that the project does not distinguish
   deuteration, so whatever an experiment or simulation with a deuterated lipid is added, we use the
   entity of a usual lipid. However, we should mention deuteration in reagent sources. Lipid synthesis
   is usually described in the paper, and it is often quite large, so we do not copy the whole
   synthesis methodology to the metadata, but we mention ``POPE: ethanolamine group is
   alpha-deuterated, synthesized according to (Vanco, 1983) from POPA (Avanti Polar Lipids)``.

   Solubles, which we explicitly add, we also mention in sources, and for water too, e.g.:

   .. code-block:: yaml

      REAGENT_SOURCES:
         TOCL: Avanti Polar Lipids
         water: Type I Milli-Q, degassed by argon bubbling
         TMACl: Sigma-Aldrich

   We do mention water as well. It can be Type I, distilled, degassed, deuterium-depleted. It’s
   important for reproducibility, and we appreciate a data curator mentioning it.

.. _metanmr_guidelines:

Guidelines to fill experiment metadata
--------------------------------------

#. **Knowledge of real temperature**

   It is not clear if the experiment is performed with knowledge of the real temperature. The pulse
   program often imposes heating inside the sample rotor, which cannot be measured directly. The
   correction should be applied, but it is sometimes unclear whether it was. This must be reflected
   in the metadata field called ``T_RF_HEATING``. If you cannot find any signs of temperature
   correction, use "UNKNOWN".

#. **Knowledge of the sign**

   99.99% of experiments don't give a sign of the order parameter, just the absolute value, and we
   should guess them when we add. The easiest way is to inherit this guess from MD. There are
   experiments where signs are measured, then we write them explicitly in ``SIGN_MEASURED`` field.

#. **Details of pulse sequence**

   Modern techniques like R-PDLF have a lot of settings in their pulse sequences. We want to reflect
   some of them in ``DETAILS`` field, but don't make it too large. If some papers are cited, we can
   cite them as well in DETAILS. Let's keep this field to a few-line paragraph.
