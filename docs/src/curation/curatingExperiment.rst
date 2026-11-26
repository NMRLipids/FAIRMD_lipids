Experiment Curation
===================

Submission structure and file naming
-------------------------------------

#. The experiment should be stored in the appropriate path.

#. The experiment should contain metadata. CI will run automatically on the pull request and check
   the fields.

#. Measurement files should have the correct names and format. NMR measurement files must use
   universal atom names.

#. If the experiment uses a molecule, the molecule should have been registered before or together
   with the experiment.

NMR metadata and measurement files
----------------------------------

It is important to verify that the experiment information is properly included, both for sample
characterisation and for the measurement method.

Characteristics of the sample are important for matching. Although many metadata fields are not
explicitly used currently by the matching algorithm, they could be useful in the special cases where
matching is performed manually. The properties are listed below in descending order of priority:

#. ~Temperature and hydration~. These characteristics can play an even more important role than accurate
   membrane composition. For example, POPC can behave more differently in 10% and 50% hydrated samples
   than in POPC:CHOL (1:3) and (1:4) samples. Check the units: hydration must be given as a water
   mass percentage.

#. ~Membrane composition~. It is of course important, however we tolerate some levels of
   uncertainty. Check that membrane composition is given in the proper units (we use molar fractions,
   not percentages).

#. ~Sample protocol~. The minimum information must be provided, but the block should not be very
   long. Pay attention to details that point to ambiguity in the sample composition. For example,
   the way of achieving a lamellar phase indicates whether the authors controlled hydration or
   simply measured the resulting state. The same applies to ionic composition, which depends strongly
   on how the sample with anionic lipids was created.

#. ~Solution composition and counterions~. Please verify that pH and ionic strength are
   reported correctly. Very often, they are unknown and it should be written in the metadata
   explicitly despite authors prefer very often to hide amiguity in the methodology section and make
   their experiments look very well controlled. Our task is completely different - ambiguous data
   should look as ambigous as possible.

#. The ~Additional molecules~ and ~Reagent sources~ sections should report details and purity.
   Additional molecules may be absent, but reagent sources should be completed where possible. Chemical
   and water purity must be mentioned. Also check that deuteration is mentioned here if 2H
   NMR was used or deuterated water was used to increase S/N ratio in 1H NMR. Pay attention to
   chelators as they are very important to ensure the absence of divalent cations, which may have
   very high affinity to phosphate headgroups.

For NMR-related fields the order to check is:

#. ~Method~, especially 2H or 1H-13C. This is extremely important to mention because it immediately
   gives a clue about where the ambiguity lies.

#. ~Temperature ambiguity~ must be mentioned. As temperature is such an important parameter, we must
   be aware of its uncertainty.

#. ~Details~ must be present and should report both pulse-sequence details (which can be referenced)
   and the assignment strategy.

#. ~Instrument~ should be named with the field strength.

Measurement files
-----------------

#. Measurement files must have valid syntax. CI should report malformed JSON files.

#. Records in one file must contain complete data for all similar C-H pairs matching the data point.
   An end CH3 group should have three records with the same number: one carbon with each of the
   protons.

#. Order parameters must have a sign, even when the sign is guessed.
