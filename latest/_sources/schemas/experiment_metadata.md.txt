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
| TOTAL_HYDRATION | Total hydration of the system (water mass %) |
| PH | pH of the system |
| PH_METHOD | Method of pH setting or measurement (buffer / measurement) |
| REAGENT_SOURCES | Description of lipid reagents. Source, purity, etc. |
| SAMPLE_PROTOCOL | Protocol for sample preparation (in free format, references are welcome). |

**NMR specific metadata**

|             key           | description                             |
|---------------------------|-----------------------------------------|
| T_RF_HEATING | Correction of the temperature according to RF heating |
| INSTRUMENT | Instrument name (including field) |
| METHOD | Method of OP measurement "abbr1:abbr2" (see detailed explaining below) |
| SIGN_MEASURED | NONE or S-DROSS |
| DETAILS | Description of NMR experiment type |

**X-ray specific metadata**

|             key           | description                             |
|---------------------------|-----------------------------------------|
| SOURCE | Source name |
| LAMBDA | Source wavelength or range |
| QRANGE | Scattering detection range (Q-range) |
| DETECTOR | Detector type |
| DISTANCE | Distance to detector (m) |
| DATATYPE | Data type (batch or SEC) |
| EXPOSURE | Exposure time per frame (s) |
| FRAMES | Number of frames |
| SAMPLE_TYPE | 'MLV', 'SUV', 'GUV', 'OS' (oriented sample), 'BIC' |

## General fields

1. **ARTICLE_DOI**  
DOI of the original publication where the experimental data originates.

2. **DATA_DOI**  
DOI of the dataset deposition with raw NMR data (e.g., nmrXive).

3. **DATA_REF**
If the dataset doesn't have DOI, we engage to add some persistent identifier or even URL if the first doesn't exist.

3. **TEMPERATURE**  
Temperature (K) of the experiment. For NMR experiment, if `NMR:T_RF_HEATING` is 'unknown' (or not given), the reported temperature from the probe is settet here. Otherwise, please insert RF-corrected temperature.

5. **MEMBRANE_COMPOSITION**  
Dictionary of molar fractions of bilayer components. For example:
```
MEMBRANE_COMPOSITION:
  POPC: 0.93
  CHOL: 0.07
```
All the molecules should be registered in the [molecular inventory](molecule_record) in the ``membrane`` subfolder.

6. **SOLUTION_COMPOSITION**  
Dictionary of solution composition of the system (mass %), main solvent is not listed:
```
SOLUTION_COMPOSITION:
  SOD: 0.5
  CLA: 0.24
  GLUCOSE: 0.1
```
All the molecules should be registered in the [molecular inventory](molecule_record) in the ``solution`` subfolder.
Do not provide whole salts! Only separated ions. Remember that the counterions of charged lipids are also part of the solution.

7. **ADDITIONAL_MOLECULES**
Dictionary of additional molecules in the format:
```
ADDITIONAL_MOLECULES:
    TFA: trifluoroacetic acid, 0.1%
    DMSO: dimethylsulfoxide, 0.1%
    DSS: sodium trimethylsilylpropanesulfonate, 0.01%
    EDTA: ethylenediaminetetraacetic acid, 0.1 mM
```
we can use INCHI-key, CAS number or just IUPAC name. If molecule is important for
the composition, it should get the metadata inside the databank and be mentioned under
`SOLUTION_COMPOSITION` instead.

8. **TOTAL_HYDRATION**  
Mass \% of water in the sample. For NMR experiment, it is better if measured by <sup>1</sup>H MAS NMR.

9. **PH**  
pH of the system (number or UNKNOWN)

10. **PH_METHOD**  
How the pH value is got: measured by pH electrode or indicator paper, measured by NMR, set by buffer.

11. **REAGENT_SOURCES**  
Which reagents are used for lipids -- should be specified for every lipid.

12. **SAMPLE_PROTOCOL**
Protocol of liposome (or OS) preparation. A description of the preparation steps and
conditions used to obtain the sample, such as lipid composition, hydration method, extrusion,
alignment procedures, and any buffers used.
For NMR sample, it is important to mention how the targeted hydration level is reached:
lyophilised powder is hydrated, liposome suspension is dehydrated, or liposome suspension
is ultracentrifugated to get lipid-rich phase.

## NMR-specific fields

All the following fields are subfields of `NMR:` block.

1. **INSTRUMENT**  
Name of the instrument and field strength.

2. **METHOD**  
A field identifying the NMR method used (string formed as METHOD:SUBMETHOD, e.g., "2H:QE").
    - Variants for METHOD: *"2H", "CDLF", "PDLF"*  
      Two main methods are <sup>2</sup>H-NMR and <sup>1</sup>H-<sup>13</sup>C SLF (separate local field)
      NMR experiments which can be either CDLF (Carbon-detected local field) or PDLF (Proton-DLF).  
    - Sub-Method for "2H": *"SP" | "QE" | "see_comments"*  
      For <sup>2</sup>H NMR, the submethod used is either "single pulse" or "quadrupolar echo".
    - Sub-Method for "CDLF": *"REDOR" | "DIPSHIFT" | "recDIPSHIFT" | "see_comments"*  
      For CDLF method, the variants could be Rotational-Echo Double-Resonance (REDOR),
      Dipolar-Coupling chemical shift correlation (DIPSHIFT), or recoupled DIPSHIFT (recDIPSHIFT).
    - Sub-Method for "PDLF": *"DROSS" | "R18_1^7" (or other numbers characterizing R-type sequence) | "see_comments"*  
      For PDLF method, subvariants could use dipolar recoupling on-axis with scaling and shape preservation (DROSS),
      or R-type recoupling (recoupling using symmetry-based pulse sequences)

3. **SIGN_MEASURED**  
Method name  (e.g. S-DROSS) if order parameter sign was measured, NONE otherwise.

4. **T_RF_HEATING**  
How RF heating is dealt (UNKNOWN / measured / guessed)

5. **DETAILS**  
Links to the pulse sequence, corresponding paper and precise parameters if important.
Obligatory explanation if **NMR:METHOD** uses "see_comments" for SUBMETHOD.

## Scattering-specific fields

All the following fields are subfields of `XRAY:` block.

1. **SOURCE**
X-ray source description. Name of the core facilities or instrument name if laboratory source was used. Name of beamline and source if synchrotron data (e.g. EMBL P12, PETRA III).

2. **LAMBDA**
Source wavelength or range. Wavelength (and/or range) of the X-ray beam used, with units (e.g., Ångstroms).

3. **QRANGE**
Scattering detection range (Q-range). The accessible scattering vector range, typically given in 1/Å.

4. **DETECTOR**
Detector type (e.g., CCD camera, PILATUS, or other detector model).

5. **DISTANCE**
Distance to detector (m). The separation between the sample and the detector, in meters.

6. **DATATYPE**
Measurement data type. Batch mode or size-exclusion chromatography (SEC) mode.

7. **EXPOSURE**
Exposure time per frame. The total data acquisition time per measurement frame, usually given in seconds.

8. **FRAMES**
Number of frames collected for dataset.

9. **SAMPLE_TYPE**
'MLV', 'SUV', 'GUV', 'OS' (oriented sample), 'BIC'. The type of sample used, with definitions: MLV (multilamellar vesicles), SUV (small unilamellar vesicles), GUV (giant unilamellar vesicles), OS (oriented sample), BIC (bicelles).

