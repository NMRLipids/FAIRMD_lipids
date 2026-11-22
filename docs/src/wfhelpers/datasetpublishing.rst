Dataset publishing workflow
===========================

The developer-only :file:`developer/gen_op_dataset.py` script generates the three
HDF5 datasets used for publication. It does not upload data. Select one or more
outputs with ``--exps``, ``--sims`` and ``--qual``; use ``--output-dir`` to keep
generated files out of the repository root:

.. code-block:: bash

   python developer/gen_op_dataset.py --exps --sims --qual --output-dir build/datasets

The generated files are ``exp-op-dataset.h5``, ``sims-op-dataset.h5`` and
``atomic-opq-dataset.h5``. Developer workflow tests live in
:file:`tests/develop/test_dataset_publishing.py` and can be run with:

.. code-block:: bash

   pytest tests/develop/test_dataset_publishing.py
