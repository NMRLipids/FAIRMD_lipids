(readmesimu)=
# Simulation metadata (README.yaml)

Each simulation in the FAIRMD Lipids is assigned with a `README.yaml` metadata file which contains all the essential information of the simulation. These files are created from the manually contributed [info.yaml](info_files) files as described in [Adding simulations](addSimulation). The `README.yaml` files are stored as described in [Data organisation](dbstructure). 

You can view examples in the [BilayerData GitHub repository](https://github.com/NMRLipids/BilayerData/tree/main/Simulations). 

README files contain information that is manually entered into [info.yaml](info_files) files and automatically exctracted information by the [fmdl_add_simulation](add_simulation_py) program. 
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
DATEOFRUNNIG | Date when added into the databank | User-given 
DIR\_WRK | Temporary local working directory | Deprecated 
UNITEDATOM\_DICT | Hydrogen information for united atom simulations | User-given
TYPEOFSYSTEM | Lipid bilayer or something else | User-given
SYSTEM | System description in the free text format | User-given
TEMPERATURE | Temperature of the simulation. | User-given 
COMPOSITION | Molecules' names, mappings, number. | Mixed*
||
PUBLICATION | Reference to a publication(s) related to the data. | User-given
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
TRJLENGTH | Lenght of the trajectory (ps). |  Autofilled
NUMBER\_OF\_ATOMS | Number of atoms in the simulation. |  Autofilled
EXPERIMENT | Potentially connected experimental data | Gen by tools
||
ID | Unique integer identifier of the simulation | Gen by CI/CD
------------------

_\* --- names and mappings are set by user whereas number of residues is autofilled_

**Fields description**

1. **DOI** (compulsory)  
Give the DOI identity for the location of simulation files. 
Current databank works only for the data in [Zenodo](https://www.zenodo.org), but other potential sources are may be implemented in the future. 
Note that the DOI must point to a specific version of dataset in Zenodo. DOIs pointing to all versions of certain dataset do not work.

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
Lenght of the trajectory (ps).

26. **TEMPERATURE**  
Temperature of the simulation.

27. **NUMBER\_OF\_ATOMS**  
Total number of atoms in the simulation.

28. **DATEOFRUNNIG**  
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
