"""
Elements assigneer functions.

Helps Molecule to define brutto formula from existing mapping files and element guesser.
"""

import json
import os
import re

import MDAnalysis as mda
import periodictable

from fairmd.lipids import FMDL_DATA_PATH
from fairmd.lipids.core import System


def uname2element(mapping_name: str) -> str:
    """Guess element name from universal atom name.

    :param mapping_name: name of the mapping, e.g. 'M_C1', 'M_G'
    :raises KeyError: if element cannot be determined from mapping name
    :return: element name, e.g. 'C', 'N', ...
    """
    # removes numbers and '_M' from the end of string
    name1 = re.sub(r"[0-9]*_M$", "", mapping_name)
    # removes M_X12Y23... sequence from the beginning
    name2 = re.sub(r"^M_([A-Z]{1,2}[0-9]{1,4})*", "", name1)

    # name2 is the atom and electrons are assigned to this atom

    if name2 == "G":  # G is a glycerol carbon so change G to C
        name2 = "C"
    elif name2 in ["X", "D"]:  # X is dummy, D is a drude particle
        return "Dummy"  # Dummy is used for unknown elements in MDAnalysis
    try:
        _ = getattr(periodictable, name2).number
    except AttributeError as e:
        msg = f"This mapping name cannot be read by our rules: {mapping_name}"
        raise KeyError(msg) from e
    else:
        return name2


def guess_elements(system: System, u: mda.Universe) -> None:
    """Assign elements to atoms in the MDAnalysis Universe based on the system's composition.

    :param system: Databank System object
    :type system: System
    :param u: MDAnalysis Universe object
    :type u: mda.Universe
    """
    for _mol, comp in system["COMPOSITION"].items():
        # get the name of molecule used in simulation files
        mol_simu_name = comp["NAME"]
        mol = system.content[_mol]

        # if lipid is split to multiple residues
        _moltype = f"moltype {comp['MOLTYPE']} and " if "MOLTYPE" in comp else ""

        try:
            ua_dict_f = os.path.join(FMDL_DATA_PATH, "lipid_json_buildH", system["UNITEDATOM_DICT"][_mol] + ".json")
            with open(ua_dict_f) as f:
                ua_dict = json.load(f)
        except (KeyError, TypeError):
            ua_dict = False

        for uname, props in mol.mapping_dict.items():
            res = str(props.get("RESIDUE", ""))
            # TODO: remove after MDAnalysis fixes the bug with "resname water" (2.11+)
            if mol_simu_name == "water":
                mol_simu_name = "water*"
            _residue = f"resname {mol_simu_name if res == '' else res} and "
            _name = f"name {props['ATOMNAME']}"
            selstr = _moltype + _residue + _name
            selection = u.select_atoms(selstr)
            if ua_dict:
                if props["ATOMNAME"] in ua_dict:
                    elname = ua_dict[props["ATOMNAME"]][0]
                    if elname == "CH":
                        elname = "CH1"
                    if elname == "CHdoublebond":
                        elname = "CH1"
                elif selection.n_atoms == 0:
                    continue  # this is an implicit atom
                else:
                    elname = uname2element(uname)
            else:
                # if not united-atom, ALL atoms must be present
                if selection.n_atoms == 0:
                    msg = f"Selection '{selstr}' did not match any atoms in the universe."
                    raise KeyError(msg)
                elname = uname2element(uname)
            selection.atoms.elements = selection.n_atoms * [elname]
        # end mapping loop
    # end molecules loop
