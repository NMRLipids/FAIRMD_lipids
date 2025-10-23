(molecule_names)=

# Universal molecule and atom names

## Molecule names
To enable automatic analyses over all simulations, universal names for molecules are defined in the NMRlipids databank as listed in the table below. These names are connected to simulation specific molecule names using the COMPOSITION dictionary in README.yaml files.

Abbreviation | Molecule name 
------------ | -------------
POPC |  1-palmitoyl-2-oleoyl-sn-glycero-3-phosphocholine
POPG |  1-palmitoyl-2-oleoyl-sn-glycero-3-phosphoglycerol
POPS | 1-palmitoyl-2-oleoyl-sn-glycero-3-phospho-L-serine
POPE | 1-palmitoyl-2-oleoyl-sn-glycero-3-phosphoethanolamine
PYPC | 1-(16:0)-2-(16:1$^\Delta9$)-sn-glycero-3-phosphocholine
PAzePCprot | 1-palmitoyl-2-azelaoyl-sn-glycero-3-phosphocholine protonated
PAzePCdeprot | 1-palmitoyl-2-azelaoyl-sn-glycero-3-phosphocholine deprotonated
DMPC | 1,2-dimyristoyl-sn-glycero-3-phosphocholine
DPPC | 1,2-dipalmitoyl-sn-glycero-3-phosphocholine
DPPE | 1,2-dipalmitoyl-sn-glycero-3-phosphoethanolamine
DPPG | 1,2-dipalmitoyl-sn-glycero-3-phospho-(1'-rac-glycerol) (sodium salt)
DEPC | 1,2-dierucoyl-sn-glycero-3-phosphocholine
DRPC | 1,2-(14:1$^\Delta9$)-sn-glycero-3-phosphocholine
DYPC | 1,2-(16:1$^\Delta9$)-sn-glycero-3-phosphocholine
DLPC | 1,2-dilauroyl-sn-glycero-3-phosphocholine
DLIPC| 1,2-dilinoleoyl-sn-glycero-3-phosphocholine
DOG  | 1,2-dioleoyl-sn-glycerol
DOPC | 1,2-dioleoyl-sn-glycero-3-phosphocholine
DOPE | 1,2-dioleoyl-sn-glycero-3-phosphoethanolamine
DDOPC| 1,2-didocosahexaenoyl-sn-glycero-3-phosphocholine
DOPS | 1,2-dioleoyl-sn-glycero-3-phospho-L-serine
DSPC | 1,2-distearoyl-sn-glycero-3-phosphocholine
DAPC | 1,2-diarachidonoyl-sn-glycero-3-phosphocholine
SLiPC | 1-(18:0)-2-(18:2 $^{\Delta9,12}$)-sn-glycero-3-phosphocholine 
DMTAP | 1,2-dimyristoyl-3-trimethylammonium-propane
GM1  | GM1 Ganglioside
SOPC | 1-stearoyl-2-oleoyl-sn-glycero-3-phosphocholine
POPI | 1-palmitoyl-2-oleoyl-sn-glycero-3-phosphoinositol 
SAPI | 1-stearoyl-2-arachidonoyl-sn-glycero-3-phosphoinositol
SAPI24 | 1-stearoyl-2-arachidonoyl-sn-glycero-3-phospho-(1'-myo-inositol-4',5'-bisphosphate)
SLPI | 1-stearoyl-2-lauroyl-sn-glycero-3-phosphoinositol
SDG | 1-stearoyl-2-docosahexaenoyl-sn-glycerol
SDPE | 1-stearoyl-2-docosahexaenoyl-sn-glycero-3-phosphoethanolamine
SM16 | N-palmitoyl-D-erythro-sphingosylphosphorylcholine
SM18 | N-stearoyl-D-erythro-sphingosylphosphorylcholine
TOCL | 1',3'-Bis[1,2-dioleoyl-sn-glycero-3-phospho]-glycerol
TLCL | tetralinoleoyl cardiolipin
CER  | N-palmitoyl-D-erythro-sphingosine
CER180 | N-stearoyl-D-erythro-sphingosine
CHOL | cholesterol 
DCHOL | 18,19-di-nor-cholesterol
DHMDMAB | dihexadecyldimethylammonium
DPPGK | 1,2-dioleoyl-sn-glycero-3-[phospho-rac-(3-lysyl(1-glycerol))] (lysyl-PG)
POT | potassium ion 
SOD | sodium ion 
CLA | chloride ion
CAL | calcium ion 
CES | caesium ion
C20 | n-eicosane
SOL | water 


## Universal atom names in mapping files
To enable automatic analyses over all simulations, universal atom names for each molecule are defined in the NMRlipids databank using the **mapping files**. In these files, universal atom names are connected to simulation specific atom names using python dictionaries stored in yaml file format. The first key in the mapping file dictionary is the universal atom name, second keys define the simulation specific atom name (`ATOMNAME`) and molecule fragment (`FRAGMENT:` head group, glycerol backbone, sn-1 or sn-2). For example, the beginning of the mapping file for CHARMM36 POPC looks like this:

     M_G1_M:
      ATOMNAME: C3
      FRAGMENT: glycerol backbone
    M_G1H1_M:
      ATOMNAME: HX
      FRAGMENT: glycerol backbone
    M_G1H2_M:
      ATOMNAME: HY
      FRAGMENT: glycerol backbone
    M_G1O1_M:
      ATOMNAME: O31
      FRAGMENT: glycerol backbone
    M_G1C2_M:
      ATOMNAME: C31
      FRAGMENT: sn-1
    M_G1C2O1_M:
      ATOMNAME: O32
      FRAGMENT: sn-1
    .
    .
    .

Universal atom names start with "M_" flag and ends with "_M" flag. In the actual naming convention between the flags, the first two characters define in which glycerol backbone chain the atoms attached (G1, G2 or G3), third character tells the atom type and fourth character tells the counting number from the glycerol backbone carbon. If there are hydrogens or other atoms attached to the main chain, those will be added to the end of the naming. More details can be found from [the original NMRlipids project post defining the mapping files](https://nmrlipids.blogspot.com/2015/03/mapping-scheme-for-lipid-atom-names-for.html). Examples already existing mapping files can be found from [the NMRlipids databank git](https://github.com/NMRLipids/Databank/tree/main/Scripts/BuildDatabank/mapping_files).

## Cardiolipin universal naming conventions

If you add a lipid of this family, you can stick to following namings:
```
                        G02O1
                    G01--G02--G03
                  /               \
               G23P1O1           G13P1O1
  G23P1O[2/3]--G23P1             G13P1--G12P1O[2/3]
               G23O1             G13O1
                 |                |
   G22C2O1       G23              G13         G12C2O1
   ||            |                |            ||
   G22C2--G22O1--G22              G12--G12O1--G12C2
  G22C3          |                |            G12C3
 G22C4           G21              G11           G12C4
 ..           G21O1                G11O1          ..
      G21C2O1==G21C2              G11C2==G11C2O1
                G21C3              G11C3
               G32C4              G11C4
                ..                 ..
```
