(readmesimu)=
# Simulation metadata (README.yaml)

Each simulation in the FAIRMD Lipids is assigned with a `README.yaml` metadata file which contains all the essential information of the simulation. These files are created from manually supplied `info.yaml` input files as described in [Adding simulations](addSimulation). Info files are not stored in BilayerData. The `README.yaml` files are stored as described in [Data organisation](dbstructure).

You can view examples in the [BilayerData GitHub repository](https://github.com/NMRLipids/BilayerData/tree/main/Simulations). 

README files contain information that is manually entered into `info.yaml` input files and automatically extracted information by the [fmdl_add_simulation](add_simulation_py) program.
Table below lists the manually entered compulsory and optional parameters, as well as automatically extracted information from simulation files.  

--------------------
key | description | type  
----|--------|------
DOI | DOI from where the raw data is found | User-given 
TRJ | Name of the trajectory file found from DOI |  User-given
TPR | Name of the topology file found from DOI |  User-given
SOFTWARE | Software used to run the simulation |  User-given
PREEQTIME | Pre-equilibrate time in nanoseconds. |  User-given
TIMELEFTOUT | Equilibration period in the uploaded trajectory. | User-given
DATEOFRUNNING | Date when added into the databank | User-given 
DIR\_WRK | Temporary local working directory | Deprecated 
UNITEDATOM\_DICT | Hydrogen information for united atom simulations | User-given
TYPEOFSYSTEM | Lipid bilayer or something else | User-given
SYSTEM | System description in the free text format | User-given
TEMPERATURE | Temperature of the simulation. | User-given 
COMPOSITION | Molecules' names, mappings, number. | Mixed*
||
PUBLICATION | Reference to a publication(s) related to the data. Legacy: enriched entries carry it as bioschema\_properties.citation | User-given
AUTHORS\_CONTACT | Name and email of the main author(s) of the data. |  User-given
BATCHID | Identifier for a series of related simulations |  User-given
SOFTWARE\_VERSION | Version of the used software |  User-given
FF | Name of the used force field |  User-given
FF\_SOURCE | Source of the force field parameters |  User-given
FF\_DATE |  Date when force field parameters were accessed |  User-given
FF{molname} | Molecule specific force field information. |  User-given
CPT | Name of the Gromacs checkpoint file. |  User-given
LOG | Name of the Gromacs log file. |  User-given
TOP | Name of the Gromacs top file. |  User-given
GRO | Name of the Gromacs gro file. |  User-given
EDR | Name of the Gromacs edr file. |  User-given
||
TRAJECTORY\_SIZE | Size of the trajectory file in bytes | Autofilled 
TRJLENGTH | Length of the trajectory (ps). |  Autofilled
NUMBER\_OF\_ATOMS | Number of atoms in the simulation. |  Autofilled
EXPERIMENT | Potentially connected experimental data | Gen by tools
||
ID | Unique integer identifier of the simulation | Gen by CI/CD
bioschema\_properties | Bioschemas Dataset description of the entry, see [below](bioschemasimu) | Gen by tools
------------------

_\* --- names and mappings are set by user whereas number of residues is autofilled_

**Fields description**

1. **DOI** (compulsory)  
Give the DOI identity for the location of simulation files. 
Current databank works only for the data in [Zenodo](https://www.zenodo.org), but other potential sources are may be implemented in the future. 
Note that the DOI must point to a specific version of dataset in Zenodo. DOIs pointing to all versions of certain dataset do not work.

`DOI` is **compulsory for every simulation entry**, and stays so: it is in the `required` list
of the README schema, and without it there is nowhere for the trajectory to be fetched from. The
metadata enrichment tooling tolerates its absence — it composes a description that simply does
not name a source rather than failing — but that is defensive handling for half-built entries,
not permission to contribute a simulation without a deposition. Entries added through
[`fmdl_add_simulation`](add_simulation_py) always have one.

2. **TRJ** (compulsory)  
Give the name of the trajectory file that is found from the DOI given above.

3. **TPR** (compulsory)  
Give the name of the file with topology information (tpr file in the case of Gromacs) that is found from the DOI given above.

3. **SOFTWARE** (compulsory)  
Give the name of software used to run the simulation. The options are GROMACS, AMBER, NAMD, CHARMM and OPENMM. So far, only simulations run with GROMACS are accepted by the script.

4. **PREEQTIME** (compulsory)  
Give the time simulated before the uploaded trajectory in nanoseconds. For example, if you upload 100-200 ns part of total 200 ns simulation, this should value should be 100.

5. **TIMELEFTOUT** (compulsory)  
Give the time that should be considered as an equilibration period in the uploaded trajectory. Frames before the give time will be discarded in the analysis. 
For example, if you upload 0-200 ns part of total 200 ns simulation where the first 100 ns should be considered as an equilibration, this value should be 100.

6. **COMPOSITION** (compulsory)  
Information about the composition (i.e., the number of molecules) of each simulation is stored in COMPOSITION as python dictionary format.
As an input, the COMPOSITION requires the universal name for each molecule present in the simulation (for definitions, see [Universal molecule and atom names](molecule_names)) as the first keys.
For each molecule, another dictionary containing the simulation specific residue name (NAME) and the name of the mapping file (MAPPING) needs to be then defined.
Mapping file defines the universal atom names for each molecule, for details see [Universal molecule and atom names](molecule_names)).
For example, the COMPOSITION dictionary input for [a system](https://doi.org/10.5281/zenodo.259392) containing POPC, Cholesterol (CHOL), water (SOL), sodium (SOD), and chloride (CLA) is given as:

```
    COMPOSITION:
     CHOL:                                      # universal molecule name
      NAME: CHL1                                # simulation specific molecule name
      MAPPING: mappingCHOLESTEROLcharmm.txt     # name of the mapping file
     POPC:
      NAME: POPC
      MAPPING: mappingPOPCcharmm.txt
     SOL:
      NAME: TIP3
      MAPPING: mappingTIP3PCHARMMgui.txt
     SOD:
      NAME: SOD
      MAPPING: mappingSOD.txt
     CLA:
      NAME: CLA
      MAPPING: mappingCLA.txt
```

When running [fmdl_add_simulation](add_simulation_py), the numbers of molecules are added in additional ``COUNT`` keys into dictionaries of each molecule, see COMPOSITION (ouput) below. 

8. **DIR\_WRK** (deprecated)  
Give the path of the working directory in your local computer. The trajectory and topology files will be downloaded to this trajectory, and temporary files created during processing will be stored here. 

9. **UNITEDATOM\_DICT** (compulsory for united atom trajectories)  
Order parameters from united atom simulations are calculated using the [buildH code](https://github.com/patrickfuchs/buildH). 
For united atom simulations, you need to tell how hydrogens are added based on definitions in 
the [JSON UA dictionary](uadic_files).
In the case of an all atom simulation, UNITEDATOM is left empty.

10. **TYPEOFSYSTEM** (compulsory)  
Tell whether the system is lipid bilayer or something else. Only lipid bilayers are currently supported, but other systems will be included in the future.

11. **SYSTEM** (compulsory)  
Give description of system in free format. For example ''POPC with cholesterol at 301K''.

12. **PUBLICATION**  
Give reference to a publication(s) related to the data.
This is a legacy field: in enriched entries the publication references are carried by
`bioschema_properties.citation` as a list of bare DOIs.

13. **AUTHORS\_CONTACT** (compulsory) 
Give the name and email of the main author(s) of the data.

14. **BATCHID**  
Give an identifier for a series of related simulations e.g. NVT simulations with the same composition but different volume carried out to determine a surface pressure-area per lipid isotherm.

15. **SOFTWARE\_VERSION**  
Give the version of the software used.

16. **FF** (compulsory)
Give the name of the force field used used in the simulation.

17. **FF\_SOURCE**  
Describe the source of the force field parameters. For example, CHARMM-GUI, link to webpage where parameters were downloaded, or citation to a paper.

18. **FF\_DATE**  
Give the date when parameters were accessed or created. The format is day/month/year.

19. **Individual force field names for molecules**
TODO: probably doesn't work currently!
In some cases special force fields are used for certain molecules. For example, non-standard parameters for ions or other molecules have been used. These can be specified giving forcefield names separately for individual molecules. These can be given as parameters named as
FF+{abbreviation from table above}, i.e., FFPOPC, FFPOT, FFSOL etc.

20. **CPT** (Gromacs)  
Give the name of the Gromacs checkpoint file that is found from the DOI given above.
CPT stands for the name of the cpt file. 

21. **LOG** (Gromacs)  
Give the name of the Gromacs log file that is found from the DOI given above.

22. **TOP** (Gromacs)  
Give the name of the Gromacs top file that is found from the DOI given above.

23. **EDR** (Gromacs)  
Give the name of the Gromacs edr file that is found from the DOI given above.

24. **TRAJECTORY\_SIZE** 
Size of the trajectory file in bytes.

25. **TRJLENGTH**  
Length of the trajectory (ps).

26. **TEMPERATURE**  
Temperature of the simulation.

27. **NUMBER\_OF\_ATOMS**  
Total number of atoms in the simulation.

28. **DATEOFRUNNING**  
Date when added into the databank.

29. **EXPERIMENT**  
Potentially connected experimental data.

30. **COMPOSITION** (output)  
When adding a simulation into FAIRMD Lipids with the [fmdl_add_simulation](add_simulation_py), numbers of molecules are automatically calculated and stored into ``COUNT`` keys for each molecule in the COMPOSITION dictionary. Number of lipids are calculated separately for both membrane leaflets. For example, the result for the simulation with COMPOSITION input exemplified above is:

```
        COMPOSITION:
         CHOL:
          NAME: CHL1
          COUNT:
           - 25
           - 25
          MAPPING: mappingCHOLESTEROLcharmm.yaml
         POPC:
          NAME: POPC
          COUNT:
           - 100
           - 100
          MAPPING: mappingPOPCcharmm.yaml
         SOL:
          NAME: TIP3
          COUNT: 9000
          MAPPING: mappingTIP3PCHARMMgui.yaml
         SOD:
          NAME: SOD
          COUNT: 21
          MAPPING: mappingSOD.yaml
         CLA:
          NAME: CLA
          COUNT: 21
          MAPPING: mappingCLA.yaml
```

31. **ID**  
Unique numeric identifier assigned to the simulation. It is generated by the CI/CD pipelines of the databank.

(bioschemasimu)=
## The `bioschema_properties` block

Enriched entries carry an extra top-level `bioschema_properties:` block, a machine-readable
description of the entry following the
[Bioschemas Dataset profile 1.0-RELEASE](https://bioschemas.org/profiles/Dataset/1.0-RELEASE),
from which schema.org JSON-LD is published. It is **written by the metadata enrichment tooling**
from the deposition record (DataCite/Zenodo) — contributors do not fill it in, and it is
optional as far as the schema is concerned.

Only properties of that Bioschemas profile are accepted (`name`, `alternateName`, `description`,
`identifier`, `keywords`, `license`, `url`, `citation`, `creator`, `datePublished`,
`distribution`, `isBasedOn`, `measurementTechnique`, `publisher`, `variableMeasured`, `sameAs`,
…), plus `dct:`-prefixed DCMI terms. Anything else is rejected, so a typo in a property name is
caught rather than silently published. One extra is tolerated for now and will be remapped:
`accessRights` (properly DCMI `dct:accessRights`).

Publication references live here as `citation`, a list of bare DOIs, superseding the legacy
top-level `PUBLICATION` field. The block also carries `_source`, local bookkeeping recording
which API the record came from and when it was retrieved; it is not a schema.org term and is
stripped before serialising JSON-LD.

Dates (`datePublished`) are `YYYY`, `YYYY-MM` or `YYYY-MM-DD` and **must stay quoted** in the
YAML — an unquoted `YYYY-MM-DD` is parsed as a date object and then fails validation as a
non-string.

`name` and `description` are **composed from the entry's own fields** — composition,
temperature, force field, engine, trajectory length — and never taken from the deposition
record. A Zenodo title names the deposition, which up to 27 entries share, so it does not tell
them apart; the free-text `SYSTEM` is kept as `alternateName`, and the deposition title as
`isPartOf.name`. Every composed title ends in the databank `ID` in brackets, which is unique.

```yaml
bioschema_properties:
  name: Molecular dynamics trajectory of a POPC bilayer at 313 K (CHARMM36, GROMACS 5.0.4,
    200 ns) [NMRlipids simulation 566]
  alternateName: 200POPC_9000SOL_313K
  description: Molecular dynamics trajectory of a lipid bilayer of 200 POPC
    (1-palmitoyl-2-oleoyl-sn-glycero-3-phosphocholine) at 313 K. Simulated with CHARMM36 in
    GROMACS 5.0.4 for 200 ns (40000 atoms). Deposited as NMRlipids Databank simulation 566 and
    available from https://doi.org/10.5281/zenodo.4040423.
  sameAs: https://doi.org/10.5281/zenodo.4040423
  datePublished: '2020-09-21'
  license:
    spdx: CC-BY-4.0
    name: Creative Commons Attribution 4.0 International
    url: https://spdx.org/licenses/CC-BY-4.0.html
    sameAs: https://creativecommons.org/licenses/by/4.0/legalcode
  publisher: Zenodo
  creator:
  - name: Ollila
    identifier: https://orcid.org/0000-0002-8135-3562
  citation:
  - 10.1021/acs.jctc.5b00935
  keywords:
  - '@type': DefinedTerm
    name: Molecular dynamics
    termCode: topic_0176
    inDefinedTermSet: http://edamontology.org
    url: http://edamontology.org/topic_0176
  - POPC
  measurementTechnique:
  - Molecular dynamics simulation (gromacs 5)
  - CHARMM36 force field
  variableMeasured:
  - area per lipid
  - C-H bond order parameter
  distribution:
  - '@type': DataDownload
    contentUrl: https://doi.org/10.5281/zenodo.4040423
    name: 200POPC_9000SOL_313K
    contentSize: 3979765856
    encodingFormat: application/x-xtc
    hasPart:
    - run.xtc
    - run.tpr
    - topol.top
  isBasedOn:
  - experiments/FormFactors/10.1016/j.bbamem.2011.07.022/11
  isPartOf:
    '@type': Dataset
    '@id': https://doi.org/10.5281/zenodo.4040423
    identifier: 10.5281/zenodo.4040423
    url: https://doi.org/10.5281/zenodo.4040423
    name: Simulation trajectories of POPC bilayers
    publisher: Zenodo
  accessRights: openAccess
  _source:
    api: datacite
    doi: 10.5281/zenodo.4040423
    spdxLicenseList: 3.28.0
    retrieved: '2026-09-08'
```
