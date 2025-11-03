#!/usr/bin/env python3

"""
A helper program to add SMILEIDX to mapping files.

This program align MD molecules from the Universe to an rdkit molecule
created by SMILES. Can be used with a care to semi-manualy curate mapping
files, SMILES, and add the alignments for only heavy atoms. Further should
be converted into a mature CI/CD script.
"""

import os
import re
import subprocess
import sys
from copy import deepcopy

import MDAnalysis as mda
import numpy as np
import yaml
from rdkit import Chem, RDLogger
from rdkit.Chem import MolStandardize

import fairmd.lipids as dlb
import fairmd.lipids.core as dlc
import fairmd.lipids.databankio as dbio
import fairmd.lipids.settings.elements as elements

lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)


def mk_universe_compact(s: dlc.System) -> mda.Universe:
    """Make universe based on just TPR and GRO records"""
    tpr_path = None
    if "TPR" in s and s["TPR"] is not None:
        tpr = s["TPR"][0][0]
        tpr_url = dbio.resolve_download_file_url(s["DOI"], tpr, validate_uri=True)
        tpr_path = os.path.join(dlb.FMDL_SIMU_PATH, s["path"], tpr)
        if not os.path.isfile(tpr_path):
            dbio.download_resource_from_uri(tpr_url, tpr_path)
    elif "PDB" in s and s["PDB"] is not None:
        msg = "PDB is not implemented yet"
        raise NotImplementedError(msg)
    else:
        msg = "System without topology."
        raise RuntimeError(msg)
    if "GRO" in s and s["GRO"] is not None:
        gro = s["GRO"][0][0]
        gro_url = dbio.resolve_download_file_url(s["DOI"], gro, validate_uri=True)
        gro_path = os.path.join(dlb.FMDL_SIMU_PATH, s["path"], gro)
        if not os.path.isfile(gro_path):
            dbio.download_resource_from_uri(gro_url, gro_path)
    elif tpr_path is not None:
        gro_path = tpr_path + ".gro"
        subprocess.run(["gmx", "editconf", "-f", tpr_path, "-o", gro_path, "-pbc"])
    else:
        msg = "Cannot generate GRO file"
        raise RuntimeError(msg)

    print(tpr_path, gro_path)
    return mda.Universe(tpr_path, gro_path)


def get_1mol_selstr(comp_name: str, mol_obj: dlc.Molecule) -> str:
    """Return selection string for a single molecule"""
    res_set = set()
    try:
        for atom in mol_obj.mapping_dict:
            res_set.add(mol_obj.mapping_dict[atom]["RESIDUE"])
    except (KeyError, TypeError):
        res_set = {comp_name}
    return "resname " + " or resname ".join(sorted(res_set))


def get_brutto_formula(eorder: str, agrp: mda.AtomGroup, charge: float = 0) -> str:
    """Get brutto formula (according to element order) of neutralized form"""
    ans = ""
    for e in eorder:
        ans += e
        n_ = (agrp.atoms.elements == e).sum()
        if e == "H" and charge < 0:
            n_ -= int(charge)
        ans += "" if n_ == 1 else str(n_)
    return ans


def compare_neutralized(a: Chem.rdchem.Mol, b: Chem.rdchem.Mol):
    """Compare neutralized forms of molecules"""
    a_ = MolStandardize.rdMolStandardize.ChargeParent(a)
    b_ = MolStandardize.rdMolStandardize.ChargeParent(b)
    aib = a_.HasSubstructMatch(b_)
    bia = b_.HasSubstructMatch(a_)
    return aib and bia


def find_uname(x, atom_name, res_name) -> str:
    for k, v in x.mapping_dict.items():
        if v["ATOMNAME"] == atom_name:
            if "RESIDUE" in v and res_name != v["RESIDUE"]:
                continue
            return k
    msg = f"Atom {atom_name} / {res_name} not found"
    raise KeyError(msg)


ss = dlc.initialize_databank()

with open("done.txt") as fd:
    done_ids = list(map(int, fd.readlines()))

# list of deeply defective systems
done_ids += [151]

print("Loading current state..")
print(done_ids)

for s in ss:
    print(s)
    if s["ID"] in done_ids:
        print(" -- Already done. Skipping.")
        continue
    if "TPR" not in s or s["TPR"] is None:
        print(" -- non-TPR currently not implemented. Skipping.")
        continue
    if s.get("UNITEDATOM_DICT", False):
        print(" -- Cannot work with UA")
        continue
    lip_names = set(s.content.keys()).intersection(dlb.core.lipids_set.names)
    print("LIPID COMPOSITION: ", lip_names)

    while lip_names:
        cur_lip = lip_names.pop()
        print("Current lipids ", cur_lip)
        mol_obj = s.content[cur_lip]
        # head dict
        print({x: mol_obj.mapping_dict[x] for x in list(mol_obj.mapping_dict.keys())[:5]})

        u = mk_universe_compact(s)

        # use internal element guesser
        u.guess_TopologyAttrs(force_guess=["elements"])
        elements.guess_elements(s, u)

        # select all molecules in the system
        sel_str = get_1mol_selstr(s["COMPOSITION"][cur_lip]["NAME"], mol_obj)
        print(sel_str)
        all_atoms_of_mol = u.select_atoms(sel_str)
        print(all_atoms_of_mol)

        # start checking the consistency
        metadata_bformula = (
            mol_obj.metadata["bioschema_properties"]["molecularFormula"].replace("+", "").replace("-", "")
        )
        eorder = re.sub(r"\d+", "", metadata_bformula)
        metadata_charge = float(mol_obj.metadata["NMRlipids"]["charge"])
        metadata_mweight = float(mol_obj.metadata["bioschema_properties"]["molecularWeight"]) + metadata_charge
        smiles = mol_obj.metadata["bioschema_properties"]["smiles"]
        mol_from_smiles = Chem.MolFromSmiles(smiles)

        molecules = all_atoms_of_mol.groupby("molnums")
        last_good_mol = None
        for _, mol_ in molecules.items():
            if u.atoms.select_atoms(f"molnum {_}") != mol_:
                continue
            try:
                mol_from_md = mol_.atoms.convert_to("rdkit")
                last_good_mol = mol_
            except Chem.AtomValenceException:
                print(f"Molecule {_} has bad conformation. Trying another one.", file=sys.stderr)
                last_good_mol = None
                continue
            print(".", end="", flush=True)
            mass_close = np.isclose(mol_.masses.sum(), metadata_mweight, atol=0.1)
            if not mass_close:
                # it can be because of incorrect average isotopic mass
                print(
                    f"Masses do not correspond to each other: {_} {mol_.masses.sum()} / {metadata_mweight}.",
                    file=sys.stderr,
                )
            cur_brutto = get_brutto_formula(eorder, mol_, metadata_charge)
            assert cur_brutto == metadata_bformula, (
                f"Brutto formulas do not correspond to each other {cur_brutto} / {metadata_bformula}"
            )
            assert compare_neutralized(mol_from_smiles, mol_from_md), "SMILES != MD"
        # check-ups are done

        # new mapping
        new_mapping_dict = deepcopy(mol_obj.mapping_dict)
        # make no-explicit-H-mol and check match-to-smiles
        mm = Chem.RemoveHs(mol_from_md)
        mtch = mm.GetSubstructMatch(mol_from_smiles)
        assert not (set(range(mol_from_smiles.GetNumHeavyAtoms())) - set(mtch)), (
            "Match string should be a permutation of [1,...,N] where N is the number of heavy atoms"
        )

        rdatit = mol_from_md.GetAtoms()
        rdatit2 = mm.GetAtoms()
        mtch_iter = iter(np.argsort(mtch))
        for mdat in last_good_mol.atoms:
            a = next(rdatit)
            if a.GetAtomicNum() == 1:
                print(a.GetAtomicNum(), mdat.name)
            else:
                b = next(rdatit2)
                match_in_smile = next(mtch_iter)
                un = find_uname(mol_obj, mdat.name, mdat.resname, s["COMPOSITION"][cur_lip]["NAME"])
                print(un, a.GetAtomicNum(), mdat.name, b.GetAtomicNum(), match_in_smile)
                new_mapping_dict[un]["SMILEIDX"] = int(match_in_smile)

        we_have_idx = False
        for k, v in new_mapping_dict.items():
            if "SMILEIDX" in v:
                if not we_have_idx and "SMILEIDX" in mol_obj.mapping_dict[k]:
                    we_have_idx = True
                if we_have_idx:
                    assert "SMILEIDX" in mol_obj.mapping_dict[k], "SMILEIDX is set, but not for this atom!"
                    assert mol_obj.mapping_dict[k]["SMILEIDX"] == v["SMILEIDX"], "SMILEIDX is wrong!!"
        if not we_have_idx:
            print("Saving mapping with SMILE indexes for the first time!")
            with open(mol_obj._mapping_fpath, "w") as fd:
                fd.write(yaml.safe_dump(new_mapping_dict, sort_keys=False, indent=1))
    with open("done.txt", "a") as fd:
        fd.write("%d\n" % s["ID"])
