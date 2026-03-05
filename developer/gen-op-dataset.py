#!/usr/bin/env python3

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from unidecode import unidecode

from fairmd.lipids.experiment import ExperimentCollection, OPExperiment

H5_MASTER = "exp-op-dataset.h5"


def store_to_hdf5(df: pd.DataFrame, e: OPExperiment, lipid_name: str):
    _mcontent = e.membrane_composition(basis="molar")
    mcontent = {"name": [], "inchikey": [], "fraction": []}
    for lname, frac in _mcontent.items():
        k = lname
        ik = e.lipids[lname].metadata["bioschema_properties"]["inChIKey"]
        mcontent["name"] += [k]
        mcontent["inchikey"] += [ik]
        mcontent["fraction"] += [frac]
    hydration = e.get_hydration()
    scontent = e.solution_composition(basis="molar")
    temperature = e["TEMPERATURE"]
    inchikey = e.lipids[lipid_name].metadata["bioschema_properties"]["inChIKey"]
    smiles = e.lipids[lipid_name].metadata["bioschema_properties"]["smiles"]
    # store all vars and df to the HDF5 table
    group = "E" + unidecode(e.exp_id).replace("-", "_").replace("/", "_").replace(".", "") + "__" + lipid_name

    with pd.HDFStore(H5_MASTER, "a") as store:
        # DataFrame table
        store.put(f"{group}/op_values", df, format="table", data_columns=True)
        # Metadata attributes
        op_storer = store.get_storer(f"{group}/op_values")
        op_storer.attrs["inchikey"] = inchikey
        op_storer.attrs["smiles"] = smiles

        # Metadata attributes
        store.put(f"{group}/sample_table", pd.DataFrame(mcontent), format="table", data_columns=True)
        sample_storer = store.get_storer(f"{group}/op_values")
        sample_storer.attrs["temperature"] = temperature
        sample_storer.attrs["hydration"] = hydration
        sample_storer.attrs["solution"] = ", ".join([f"{k:<25} {v * 100:>6.1f}%" for k, v in sorted(scontent.items())])


def by_c_uname(opdict: dict, uname: str) -> dict:
    """Return only OP vals of named heavy atoms. Sorted by val."""
    res = []
    for k, v in opdict.items():
        if k.split()[0] == uname:
            res.append(v)
    return sorted(res, key=lambda x: x[0])


def main() -> None:
    print("Generating OP datasets.")
    ee = ExperimentCollection.load_from_data("OPExperiment")
    for e in ee:
        print(e)
        for lname, opdict in e.data.items():
            mol = e.lipids[lname]
            smiles = mol.metadata["bioschema_properties"]["smiles"]
            rdmol = Chem.MolFromSmiles("CCO")  # Ethanol
            n_heavy_atoms = rdMolDescriptors.CalcNumHeavyAtoms(rdmol)
            print(smiles)
            smi2uname = {}
            for uname, aprops in mol.mapping_dict.items():
                if "SMILEIDX" in aprops:
                    smid = int(aprops["SMILEIDX"])
                    smi2uname[smid] = uname
            if not smi2uname:
                # NO SMILEIDX. Cannot store.
                continue
            smi2uname = dict(sorted(smi2uname.items()))
            df = pd.DataFrame(
                {
                    "id": np.full(len(opdict), -1, dtype=np.int64),  # int64, init NaN
                    "val": np.full(len(opdict), np.nan),
                    "err": np.full(len(opdict), np.nan),
                }
            )
            cur_row = 0
            for id_, uname in smi2uname.items():
                for opdict_row in by_c_uname(opdict, uname):
                    if cur_row >= len(df):
                        break
                    df.at[cur_row, "id"] = id_
                    df.at[cur_row, "val"] = opdict_row[0]
                    df.at[cur_row, "err"] = 0.02 if len(opdict_row) == 1 else opdict_row[1]
                    cur_row += 1
            # now, we are ready to store to HDF5
            store_to_hdf5(df, e, lname)


if __name__ == "__main__":
    main()
