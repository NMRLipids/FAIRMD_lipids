"""
The simulation and experiment Bioschemas autocomplete, split by concern.

``autocomplete_expsim_metadata.py`` beside this package is the command line; the
work is here, in one module per layer, each importing only from the ones above
it:

- :mod:`constants`    -- endpoints, licences, block layout; constants.yaml
- :mod:`helpers`      -- fetching one JSON document, normalising a string
- :mod:`licenses`     -- resolving a licence statement onto the SPDX list
- :mod:`registries`   -- DataCite and CrossRef answers as Bioschemas blocks
- :mod:`fields`       -- what a ``README.yaml`` says about itself
- :mod:`inheritance`  -- what a record inherits from its kind's template
- :mod:`descriptions` -- the composed ``name`` and ``description``
- :mod:`bioschema`    -- assembling the block for one record
- :mod:`records`      -- reading and rewriting ``README.yaml``
- :mod:`checks`       -- validating what a run wrote
"""
