.. _quality

Simulation quality
==================

Quality evaluation in the FairMD databank
----------------------------------------

The simulations contained in the FairMD databank are quality evaluated against experiments.
A simulation is connected to an experimental data set if the molar concentrations of all
molecules are within :math:`\pm 3` percentage units, charged lipids have the same
counterions, and temperatures are within :math:`\pm 2\,\mathrm{K}`.
For molar concentrations of water, the exact hydration level is considered only for
systems with molar water-to-lipid ratio below 25; otherwise the systems are considered
fully hydrated.

Two types of experimental data are used in quality evaluation: SAXS form factors and
C–H order parameters. For the latter, a quality is calculated both for the entire
molecule and for molecular fragments, as described below.

The qualities for each simulation record are saved in the files described in
:ref:`simulation-quality-files`.

The simulations contained in the databank can then be ranked based on their qualities.
This is exemplified in the ranking files provided in :ref:`ranking_files`.


Quality evaluation of C–H bond order parameters
-----------------------------------------------

The quality of each C–H bond order parameter is estimated by calculating the probability
for a simulated value to lie within the error bars of the experimental value.
Because conformational ensembles of individual lipids are assumed to be independent in
a fluid lipid bilayer,

.. math::

   \frac{S_\mathrm{CH} - \mu}{s / \sqrt{n}}

has a Student’s *t*-distribution with :math:`n - 1` degrees of freedom, where
:math:`\mu` represents the real mean of the order parameter.
The probability for an order parameter from simulation to lie within experimental error
bars can be estimated as

.. math::

   P = f\!\left(\frac{S_\mathrm{CH} - (S_\mathrm{exp} + \Delta S_\mathrm{exp})}{s / \sqrt{n}}\right)
     - f\!\left(\frac{S_\mathrm{CH} - (S_\mathrm{exp} - \Delta S_\mathrm{exp})}{s / \sqrt{n}}\right)

where :math:`f(t)` is the first-order Student’s *t*-distribution,
:math:`n` is the number of independent sample points for each C–H bond (the number of
lipids in a simulation),
:math:`S_\mathrm{CH}` is the sample mean order parameter,
:math:`s` is the variance of :math:`S_\mathrm{CH}` calculated over individual lipids,
:math:`S_\mathrm{exp}` is the experimental value, and
:math:`\Delta S_\mathrm{exp}` its error.
An error of :math:`\Delta S_\mathrm{exp} = 0.02` is assumed for all experimental order
parameters.

We use the first-order Student’s *t*-distribution, which has slightly higher
probabilities for values far from the mean, to avoid numerical instabilities stemming
from significant disagreement between simulation and experimental values.
Furthermore, some force fields exhibit too slow dynamics, leading to large error bars in
the :math:`S_\mathrm{CH}` values.
Such artificially slow dynamics widen the Student’s *t*-distribution, consequently
increasing the probability of finding the simulated value within experimental error
bars.
Therefore, :math:`S_\mathrm{CH}` values with simulation error bars above the experimental
error of 0.02 are not included in the quality evaluation.

We define the average qualities for different fragments
(frag = *sn-1*, *sn-2*, *headgroup*, or *total*, with the last referring to all order
parameters within a molecule) within each lipid type in a simulation as

.. math::

   P^\mathrm{frag}[\mathrm{lipid}]
   =
   \langle P[\mathrm{lipid}] \rangle_\mathrm{frag}
   \, F_\mathrm{frag}[\mathrm{lipid}],

where :math:`\langle P[\mathrm{lipid}] \rangle_\mathrm{frag}` is the average of the
individual :math:`S_\mathrm{CH}` qualities within the fragment, and
:math:`F_\mathrm{frag}[\mathrm{lipid}]` is the percentage of order parameters for which
the quality is available within the fragment.

The overall quality of different fragments in a simulation
(frag = *tails*, *headgroup*, or *total*) is then defined as a molar-fraction–weighted
average over different lipid components:

.. math::

   P^\mathrm{frag}
   =
   \sum_\mathrm{lipid}
   \chi_\mathrm{lipid}
   P^\mathrm{frag}[\mathrm{lipid}],

where :math:`\chi_\mathrm{lipid}` is the molar fraction of a lipid in the bilayer and
*tails* refers to the average over all acyl chains.


Quality evaluation of X-ray scattering form factors
---------------------------------------------------

We use only the location of the first form factor minimum for quality evaluation.
This choice is motivated by the simulation-size dependence of the relative form factor
lobe heights and by the difficulty of automatically detecting the locations of further
minima in some experimental data sets due to noise.

First, fluctuations are filtered from the form factor data using a Savitzky–Golay filter
(window length 30 and polynomial order 1).
The first minimum at :math:`q > 0.1\,\mathrm{Å}^{-1}` is then located for both simulation
(:math:`FF_\mathrm{min}^\mathrm{sim}`) and experiment
(:math:`FF_\mathrm{min}^\mathrm{exp}`).

The quality of a form factor is defined as the Euclidean distance between the minima
locations:

.. math::

   FF_q
   =
   \left| FF_\mathrm{min}^\mathrm{sim} - FF_\mathrm{min}^\mathrm{exp} \right| \times 100.


.. _Ranking files:
   https://github.com/NMRLipids/BilayerData/blob/view-new-rankings/Ranking/
