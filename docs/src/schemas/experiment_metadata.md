(readmeexp)=
# Experiment metadata

**Global metadata table**

|             key           | description |
|---------------------------|-----------------------------------------|
| ARTICLE_DOI | DOI of the of the original publication of the experimental data. Becomes the citation when no DATA_DOI is given |
| DATA_DOI | DOI of the dataset deposition with raw NMR data. When present, this is the DOI that becomes the citation |
| DATA_REF | Reference to deposited dataset, for data without a DOI |
| DATE | Date when the data was recorded or published |
| TEMPERATURE | Temperature (K) of the experiment |
| MEMBRANE_COMPOSITION | Dictionary of molar fractions of membrane phase |
| SOLUTION_COMPOSITION | Dictionary of ion concentrations in the system |
| ADDITIONAL_MOLECULES| Dictionary of molecules which are not in the databank |
| TOTAL_HYDRATION | Total hydration of the system (water mass %) |
| PH | pH of the system |
| PH_METHOD | Method of pH setting or measurement (buffer / measurement) |
| REAGENT_SOURCES | Description of lipid reagents. Source, purity, etc. |
| SAMPLE_PROTOCOL | Protocol for sample preparation (in free format, references are welcome). |
| bioschema_properties | Generated Bioschemas Dataset block, see [below](bioschemaexp) |

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

At least one of `ARTICLE_DOI` and `DATA_DOI` has to be given, and the two decide what the
generated `bioschema_properties.citation` holds: **if a `DATA_DOI` is present, that is the DOI
that becomes the citation** — the data is what should be cited — and only when `ARTICLE_DOI` is
the sole DOI does the article become the citation.

Both are given as bare DOIs — `10.1016/j.bbamem.2011.07.022`, not `doi:10.1016/...` and not
`https://doi.org/10.1016/...`. Data that has no DOI at all is referenced with `DATA_REF` instead.

1. **ARTICLE_DOI**  
DOI of the original publication where the experimental data originates.

2. **DATA_DOI**  
DOI of the dataset deposition with raw NMR data (e.g., nmrXive).

3. **DATA_REF**
If the dataset doesn't have DOI, we engage to add some persistent identifier or even URL if the first doesn't exist.

4. **DATE**
Date in the standard format YYYY-MM-DD (e.g., 2023-08-24). A lot of date values have been automatically synchronized from the paper dates. If the data wasn't published, the date of recording should be used.

5. **TEMPERATURE**  
Temperature (K) of the experiment, so strictly positive. For NMR experiment, if `NMR:T_RF_HEATING` is 'unknown' (or not given), the reported temperature from the probe is settet here. Otherwise, please insert RF-corrected temperature.

6. **MEMBRANE_COMPOSITION**  
Dictionary of molar fractions of bilayer components, each within (0, 1]. For example:
```
MEMBRANE_COMPOSITION:
  POPC: 0.93
  CHOL: 0.07
```
All the molecules should be registered in the [molecular inventory](molecule_record) in the ``membrane`` subfolder.

7. **SOLUTION_COMPOSITION**  
Dictionary of solution composition of the system (mass %, **not mM!**), main solvent is not listed:
```
SOLUTION_COMPOSITION:
  SOD: 0.5
  CLA: 0.24
  GLUCOSE: 0.1
```
All the molecules should be registered in the [molecular inventory](molecule_record) in the ``solution`` subfolder.
Do not provide whole salts! Only separated ions. Remember that the counterions of charged lipids are also part of the solution.

8. **ADDITIONAL_MOLECULES**
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

9. **TOTAL_HYDRATION**  
Mass \% of water in the sample, so within (0, 100]. For NMR experiment, it is better if measured by <sup>1</sup>H MAS NMR.

10. **PH**  
pH of the system (number or UNKNOWN)

11. **PH_METHOD**  
How the pH value is got: measured by pH electrode or indicator paper, measured by NMR, set by buffer.

12. **REAGENT_SOURCES**  
Which reagents are used for lipids -- should be specified for every lipid.

13. **SAMPLE_PROTOCOL**
Protocol of liposome (or OS) preparation. A description of the preparation steps and
conditions used to obtain the sample, such as lipid composition, hydration method, extrusion,
alignment procedures, and any buffers used.
For NMR sample, it is important to mention how the targeted hydration level is reached:
lyophilised powder is hydrated, liposome suspension is dehydrated, or liposome suspension
is ultracentrifugated to get lipid-rich phase.

(bioschemaexp)=
## The `bioschema_properties` block

Enriched entries carry an extra top-level `bioschema_properties:` block, a machine-readable
description of the entry following the
[Bioschemas Dataset profile 1.0-RELEASE](https://bioschemas.org/profiles/Dataset/1.0-RELEASE),
from which schema.org JSON-LD is published. It is **written by the metadata enrichment tooling**
from the deposition record (CrossRef or DataCite) — contributors do not fill it in, and it is
optional as far as the schema is concerned.

Only properties of that Bioschemas profile are accepted (`name`, `description`, `identifier`,
`keywords`, `license`, `url`, `citation`, `creator`, `datePublished`, `distribution`,
`isBasedOn`, `measurementTechnique`, `publisher`, `variableMeasured`, `isPartOf`, `sameAs`, …),
plus `dct:`-prefixed DCMI terms. Anything else is rejected, so a typo in a property name is
caught rather than silently published. Two extras are tolerated for now and will be remapped:
`accessRights` (properly DCMI `dct:accessRights`) and `articleLicense` (belongs on
`isPartOf.license`).

The block also carries `_source`, local bookkeeping recording which API the record came from
and when it was retrieved. It is not a schema.org term and is stripped before serialising
JSON-LD.

Dates (`datePublished`) are `YYYY`, `YYYY-MM` or `YYYY-MM-DD` and **must stay quoted** in the
YAML — an unquoted `YYYY-MM-DD` is parsed as a date object and then fails validation as a
non-string.

```yaml
bioschema_properties:
  name: Fluid phase lipid areas and bilayer thicknesses of commonly used phosphatidylcholines
  description: Experimental X-ray scattering form factor for a lipid bilayer containing POPC ...
  sameAs: https://doi.org/10.1016/j.bbamem.2011.07.022
  datePublished: '2011-11'
  license:
    spdx: CC-BY-4.0
    name: Creative Commons Attribution 4.0 International
    url: https://spdx.org/licenses/CC-BY-4.0.html
    sameAs: https://creativecommons.org/licenses/by/4.0/
  articleLicense:
    url: http://www.elsevier.com/open-access/userlicense/1.0/
    spdx: null
  publisher: Elsevier BV
  creator:
  - name: Norbert Kucerka
  citation:
  - 10.1016/j.bbamem.2011.07.022
  keywords:
  - '@type': DefinedTerm
    name: X-ray diffraction
    termCode: topic_2828
    inDefinedTermSet: http://edamontology.org
    url: http://edamontology.org/topic_2828
  - POPC
  measurementTechnique:
  - '@type': DefinedTerm
    name: small-angle X-ray scattering
    termCode: CHMO_0000204
    inDefinedTermSet: http://purl.obolibrary.org/obo/chmo.owl
    url: http://purl.obolibrary.org/obo/CHMO_0000204
  - X-ray scattering (SUV)
  variableMeasured:
  - '@type': PropertyValue
    name: X-ray scattering form factor
    unitText: A^-1
  distribution:
  - '@type': DataDownload
    name: POPC_ULV_20Cin0D_FormFactor.json
    encodingFormat: application/json
  isPartOf:
    '@type': ScholarlyArticle
    '@id': https://doi.org/10.1016/j.bbamem.2011.07.022
    identifier: 10.1016/j.bbamem.2011.07.022
    url: https://doi.org/10.1016/j.bbamem.2011.07.022
    name: Fluid phase lipid areas and bilayer thicknesses of commonly used phosphatidylcholines
    isPartOf:
      '@type': Periodical
      name: Biochimica et Biophysica Acta (BBA) - Biomembranes
  _source:
    api: crossref
    doi: 10.1016/j.bbamem.2011.07.022
    retrieved: '2026-09-08'
```

## NMR-specific fields

All the following fields are subfields of `NMR:` block. **INSTRUMENT**, **METHOD**,
**SIGN_MEASURED** and **T_RF_HEATING** are all required whenever the block is given;
**DETAILS** is required on top of those when **METHOD** uses `see_comments`.

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

All the following fields are subfields of `XRAY:` block. **SOURCE**, **LAMBDA** and
**SAMPLE_TYPE** are required whenever the block is given.

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

