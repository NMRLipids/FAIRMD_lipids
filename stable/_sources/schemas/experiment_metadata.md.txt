(readmeexp)=
# Experiment metadata

|key | description |
| ----| -----|
|DOI | DOI of the of the original publication of the experimental data|
TEMPERATURE | Temperature of the experiment
MOLAR_FRACTIONS | Dictionary of molar fractions of bilayer components
ION_CONCENTRATIONS | Dictionary of ion concentrations in the system
TOTAL_LIPID_CONCENTRATION | Total concentration of lipid components
COUNTER_IONS | Type of counter ions if present

**Field descriptions**

1. **DOI**  
DOI of the original publication where the experimental data originates.

2. **TEMPERATURE**  
Temperature of the experiment.

3. **MOLAR_FRACTIONS**  
Dictionary of molar fractions of bilayer components. For example:
```
MOLAR_FRACTIONS:
  POPC: 0.93
  CHOL: 0.07
```

4. **ION_CONCENTRATIONS**  
Dictionary of ion concentrations of the system ([ion]= N<sub>ion</sub>*[water]/N<sub>water</sub>, where [water] = 55.5 M). For example:
```
ION_CONCENTRATIONS:
  POT: 0
  SOD: 0.1
  CLA: 0.1
  CAL: 0
```

5. **TOTAL_LIPID_CONCENTRATION**  
Total concentration of lipid components ([lipid]= N<sub>lipid</sub>*[water]/N<sub>water</sub>, where [water] = 55.5 M). If exact concentration is not known, but experiments are performed in excess water, ’full hydration’ can be given.

6. **COUNTER_IONS**  
Type of counter ions if present.
