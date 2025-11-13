(readmeexp)=
# Experiment metadata

**Global metadata table**

|             key           | description |
|---------------------------|-----------------------------------------|
| ARTICLE_DOI | DOI of the of the original publication of the experimental data|
| DATA_DOI | DOI of the dataset deposition with raw NMR data |
| DATA_REF | Reference to deposited dataset |
| TEMPERATURE | Temperature (K) of the experiment |
| MEMBRANE_COMPOSITION | Dictionary of molar fractions of membrane phase |
| SOLUTION_COMPOSITION | Dictionary of ion concentrations in the system |
| ADDITIONAL_MOLECULES| Dictionary of molecules which are not in the databank |
| PH | pH of the system |
| PH_METHOD | Method of pH setting or measurement (buffer / measurement) |
| REAGENT_SOURCES | Description of lipid reagents. Source, purity, etc. |

**NMR specific metadata**

|             key           | description                             |
|---------------------------|-----------------------------------------|
| T_RF_HEATING | Correction of the temperature according to RF heating |
| TOTAL_HYDRATION | Total hydration of the system (water mass %) |
| HYDRATION_METHOD | Method of preparing lamellar phase |
| NMR_INSTRUMENT | Instrument name (including field) |
| NMR_EXPERIMENT | Description of NMR experiment type |

**Field descriptions**

1. **ARTICLE_DOI**  
DOI of the original publication where the experimental data originates.

2. **DATA_DOI**  
DOI of the dataset deposition with raw NMR data (e.g., nmrXive).

3. **DATA_REF**
If the dataset doesn't have DOI, we engage to add some persistent identifier or even URL if the first doesn't exist.

3. **TEMPERATURE**  
Temperature (K) of the experiment. If `T_RF_HEATING` is 'no_inromation', the reported temperature from the probe is given. Otherwise, please insert RF-corrected temperature.

4. **T_RF_HEATING**  
How RF heating is dealt (no_information / measured / guessed)

5. **MEMBRANE_COMPOSITION**  
Dictionary of molar fractions of bilayer components. For example:
```
MOLAR_FRACTIONS:
  POPC: 0.93
  CHOL: 0.07
```

6. **SOLUTION_COMPOSITION**  
Dictionary of solution composition of the system (mass %):
```
ION_CONCENTRATIONS:
  SOD: 0.5
  CLA: 0.24
  GLUCOSE: 0.1
```
Do not provide whole salts! Only separated ions. Remember that the counterions of charged lipids are also part of the solution.

7. **ADDITIONAL_MOLECULES**
Dictionary of additional molecules in the format:
```
ADDITIONAL_MOLECULES:
    TFA: trifluoroacetic acid
    DMSO: dimethylsulfoxide
    DMS: trimethylsilylpropanesulfonate
```
we can use INCHI-key, CAS number or just IUPAC name. If molecule is important for the composition, it will get the metadata inside the databank.

8. **PH**  
pH of the system (number or UNKNOWN)

9. **PH_METHOD**  
How the pH value is got: measured by pH electrode or indicator paper, measured by NMR, set by buffer.

10. **REAGENT_SOURCES**  
Which reagents are used for lipids -- should be specified for every lipid.

11. **TOTAL_HYDRATION**  
Mass \% of water. Better if it is measured by 1H MAS NMR.

12. **HYDRATION_METHOD**  
Way how the targeted hydration level is reached: lyophilised powder is hydrated, liposome suspension is dehydrated, or liposome suspension is ultracentrifugated to get lipid-rich phase.

13. **NMR_INSTRUMENT**
Name of the instrument and field strength.

14. **NMR_EXPERIMENT**
Type of NMR experiment used. Specify whatever is required. Links to the pulse sequence, corresponding paper and precise parameters would be nice.
