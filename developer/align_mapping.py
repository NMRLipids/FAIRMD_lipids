#!/usr/bin/env python3
"""
A helper program to add SMILEIDX to mapping files.

This program align MD molecules from the Universe to an rdkit molecule
created by SMILES. Can be used with a care to semi-manualy curate mapping
files, SMILES, and add the alignments for only heavy atoms. Further should
be converted into a mature CI/CD script.
"""

import argparse
import os
import re
import sys
from copy import deepcopy

import MDAnalysis as mda
import numpy as np
import yaml
from rdkit import Chem, RDLogger
from rdkit.Chem import MolStandardize

import fairmd.lipids.api as dlapi
import fairmd.lipids.core as dlc
import fairmd.lipids.molecules as dlm
from fairmd.lipids.auxiliary import elements

lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)


def get_1mol_selstr(comp_name: str, mol_obj: dlm.Molecule) -> str:
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


def compare_neutralized(a: Chem.rdchem.Mol, b: Chem.rdchem.Mol) -> bool:
    """Compare neutralized forms of molecules"""
    a_ = MolStandardize.rdMolStandardize.ChargeParent(a)
    b_ = MolStandardize.rdMolStandardize.ChargeParent(b)
    aib = a_.HasSubstructMatch(b_)
    bia = b_.HasSubstructMatch(a_)
    return aib and bia


DONE_LOG_FNAME = "align-mapping-file.log"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--skip-ids",
        type=lambda s: [int(x) for x in s.split(",")],
        default=[],
        help="Comma-separated list of IDs to skip",
    )

    args = parser.parse_args()

    ss = dlc.initialize_databank()

    done_ids = []
    if os.path.isfile(DONE_LOG_FNAME):
        with open(DONE_LOG_FNAME) as fd:
            done_ids = list(map(int, fd.readlines()))
    done_ids.extend(args.skip_ids)  # add defective systems

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
        print("LIPID COMPOSITION: ", s.lipids.keys())

        u: mda.Universe | None = None

        for cur_lip, mol_obj in s.lipids.items():
            print("Current lipids ", cur_lip)
            smileids = [x["SMILEIDX"] for x in mol_obj.mapping_dict.values() if "SMILEIDX" in x]
            # head dict
            if smileids:
                print(f" -- SMILEIDX already present ({len(smileids)}). Skipping.")
                continue
            print("START WORKING -- ", s)
            if u is None:
                uc = dlapi.UniverseConstructor(s)
                uc.download_mddata(skip_traj=True)
                u = uc.build_universe()

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
                except (Chem.AtomValenceException, KeyError):
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
                    un = mol_obj.md2uan(mdat.name, mdat.resname)
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
        with open(DONE_LOG_FNAME, "a") as fd:
            fd.write("%d\n" % s["ID"])


if __name__ == "__main__":
    main()
