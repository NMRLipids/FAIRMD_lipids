.. _addnewmol:

Adding new molecule
===================

First, please do not hesitate to ask assistance via
`BilayerData GitHub issues <https://github.com/NMRlipids/BilayerData/issues>`_.

If you are adding data into the databank and your molecule(s) ``YOURMOL`` do not exist,
you need to create a new one by adding a folder :file:`membrane/YOURMOL/` inside
:py:data:`FMDL_MOL_PATH` folder. When contributing to the standard NMRlipids Data
repository, `BilayerData <https://github.com/NMRlipids/BilayerData>`_, you are adding
into :file:`Molecules/membrane` folder.

The folder should contain at least two files: metadata file :file:`metadata.yaml` and
mapping file :file:`yourmol-forcefieldname-mapping.yaml` (this name is arbitrary, but
you should make it clear).

.. code-block:: console

   ./ ├── Simulations/
      ├── experiments/
      ├── info_files/
      ├── lipid_json_buildH/
      ├── Ranking/
      └── Molecules/
          └── membrane/
              └── YOURMOL/
                  ├── yourmol-forcefield2-mapping.yaml
                  └── metadata.yaml

Metadata file creation
~~~~~~~~~~~~~~~~~~~~~~

Please provide metadata about the molecule in a :file:`metadata.yaml` in the same
subfolder of :file:`Molecules/membrane/YOURMOLECULE`. You can find the template in
:file:`SchemaValidation/Schema/metadata_template.yaml` inside the package directory.

The recommended workflow is to start from the `InChI
<https://en.wikipedia.org/wiki/International_Chemical_Identifier>`_. You can obtain the
InChI and InChIKey for your molecule via different methods. One possiblity is to use a
PDB snapshot of your trajectory (extracted by GROMACS, VMD, MDAnalysis ...) and to add
connectivty using `RDkit <https://www.rdkit.org/>`_ or `gromologist
<https://gitlab.com/KomBioMol/gromologist>`_ and convert a single molecule to InChI via
`RDkit <https://www.rdkit.org/>`_ or `OpenBabel <https://openbabel.org/index.html>`_. We
use the neutral form of the molecule. Please indicate the charge of the molecule in your
:file:`metadata.yaml` under ``charge:``.

Please ensure that the InChI represents your molecule. The named tools might deliver
imperfect results. Provided your molecule is not fully synthetic or novel, you can
obtain most of the other metadata from the `PubChem API
<https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest>`_ and the `UniChem API
<https://www.ebi.ac.uk/unichem/api/docs>`_ using the InChIKey. Mapping especially to
LIPID MAPS and SwissLipids might be incomplete in some cases and can be completed
manually.

Mapping file creation
~~~~~~~~~~~~~~~~~~~~~

The easiest way is to take similar already existing mapping file and modify that. If
atoms in a lipid belong to different residues, which is typical situation in Amber force
fields, for example see `here
<https://github.com/NMRLipids/BilayerData/blob/main/Molecules/membrane/POPC/mappingPOPClipid17.yaml>`_,
add the residue name to ``RESIDUE`` key of each atom in the mapping file. In this case,
give the name of the head group residue in the ``COMPOSITION`` dictionary in
:ref:`the README.yaml file <readmesimu>`.

The mapping file should contain all the atoms of the molecules.

Specific cases:
- TODO: mapping for UA
