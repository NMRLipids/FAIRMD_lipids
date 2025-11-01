.. _unittesting:

Project testing
===============

Testing strategy
----------------

Test function requiring ``Simulation.1`` and ``Simulation.2`` path plugged in as
:var:`fairmd.lipids.NMLDB_SIMU_PATH` must be groupped in different test_* files. Enviromental variable
selection is performed in ``conftest.py::header_module_scope`` fixture at the
test-module level.

Tests are organized via `tox <https://tox.wiki/>`_.

.. code-block:: bash

 tox -e tests-all -- tests/test_load.py

During setting up the environment, tox will replicate ``ToyData`` from ``src/data/ToyData``
to ``tests`` folder. You can remove ``tests/ToyData`` if you don't need it for debugging.

Simulation.1 and Simulation.2
-----------------------------

In ``ToyData/Simulation.1`` we collected only ``README.yaml`` files for simulations, so that testing
functions can download everything, perform analyzis over the trajectory and write files. This
behavior corresponds to the situation when *de novo* analysis is run after adding a new system.

In ``ToyData/Simulation.2``, in contrary, we put completed Databank analysis files and we will
test only those API part on this toy-Databank, which doesn't require MDAnalysis objects
to create. Some function requires only JSON files to exist -- for them we prepared
this folder.

