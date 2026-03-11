(molecule_names)=

# Universal molecule and atom names

## Molecule names
To enable automatic analyses over all simulations, universal names for molecules are defined
in the FAIRMD Lipids as listed in the table below. These names are connected to simulation
specific molecule names using the `MEMBRANE_COMPOSITION` dictionary in `README.yaml` files.
Note that the content of the table is loaded dynamically from [FAIRMD Lipids web portal](https://databank2.nmrlipids.fi).

<iframe src="https://databank2.nmrlipids.fi/lipids?items_per_page=all&embed=1" width="100%" height="400px" frameborder="0"></iframe>

## Universal atom names in mapping files
To enable automatic analyses over all simulations, universal atom names for each molecule are defined in the FAIRMD Lipids using the **mapping files**. In these files, universal atom names are connected to simulation specific atom names using python dictionaries stored in yaml file format. The first key in the mapping file dictionary is the universal atom name, second keys define the simulation specific atom name (`ATOMNAME`) and molecule fragment (`FRAGMENT:` head group, glycerol backbone, sn-1 or sn-2). For example, the beginning of the mapping file for CHARMM36 POPC looks like this:

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

Universal atom names start with "M_" flag and ends with "_M" flag. 

## Glycerolipids naming convention

In the actual naming convention between the flags, the first two characters define in which glycerol backbone chain the atoms attached (G1, G2 or G3), third character tells the atom type and fourth character tells the counting number from the glycerol backbone carbon. If there are hydrogens or other atoms attached to the main chain, those will be added to the end of the naming. More details can be found from [the original NMRlipids project post defining the mapping files](https://nmrlipids.blogspot.com/2015/03/mapping-scheme-for-lipid-atom-names-for.html). Examples already existing mapping files can be found in [the Toy Databank](https://github.com/NMRLipids/FAIRMD_lipids/tree/main/src/fairmd/lipids/data/ToyData/Molecules).

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

## Sphingolipids naming convention

Shpingosine chain of sphingolipids is named according to the numeration of sphingosine, i.e., `M_C1_M`, `M_C2_M`, etc.
Hydroxy-group is called `M_C3O1_M`, `M_C3O1H1_M`, accordingly.
N-atom is called `M_N1_M`. N-attached fatty acid is named as `M_N1C1_M`, `M_N1C2_M`, etc.
Choline/ethanolamine headgroups should be named according to glycerolipid convention.

## Cholesterol

Cholesterol and its derivatives are named according to IUPAC numbering: `M_Ci_M` and `M_CiHj_M`.
Hydroxy-group is named according to the attached carbon: `M_C3O1_M`.