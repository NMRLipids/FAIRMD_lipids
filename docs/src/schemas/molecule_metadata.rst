.. _molmetadata_schema:

Molecule metadata schema
==========================

Project stores two types of molecules: solubles and membrane molecules (which are mostly lipids).

The metadata of solubles is intentionally small. Here is the inline template:

.. code-block:: yaml

   NMRlipids:
     id: SOL
     charge: 0
     mass: 18.0

The soluble metadata is validated with :file:`soluble_schema.json`, while membrane
molecule metadata is validated with :file:`molecule_schema.json`.

The metadata of membrane molecules is rich and described in the template
:file:`schema_validation/schema/molecule_template.yaml` which is distributed with the package.

.. literalinclude:: ../../../src/fairmd/lipids/schema_validation/schema/molecule_template.yaml
   :language: yaml
