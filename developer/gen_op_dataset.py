#!/usr/bin/env python3

import argparse
import logging
import sys

import _bootstrap
from datasets.pureop import gen_op_from_exps, gen_op_from_sims, gen_opq_ds

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="Dataset Generator",
        description="""
CI helper script for delivering dataset for Kaggle in HDF5 format.
Each dataset is composed of records representing one ssNMR experiment or simulation for one lipid,
i.e. if there are two lipids in one simulation, it will be stored as two different records.

Two DataFrames are stored for each record:
1. SMILES-aligned values for each atom of the molecule
2. Composition table of membrane part of the system""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--exps", action="store_true", help="Generate DS from experiments")
    parser.add_argument("--sims", action="store_true", help="Generate DS from simulations")
    parser.add_argument("--qual", action="store_true", help="Generate DS with quality from paired simulations")
    args = parser.parse_args()

    if not args.exps and not args.sims:
        parser.print_usage()
        sys.exit(1)

    lg = logging.getLogger("cli")
    lg.setLevel(logging.INFO)
    h_stderr = logging.StreamHandler(sys.stderr)
    h_stderr.setLevel(logging.INFO)
    lg.addHandler(h_stderr)

    if args.exps:
        e_ok, e_fail = gen_op_from_exps(lg)
    if args.sims:
        s_ok, s_fail = gen_op_from_sims(lg)
    if args.qual:
        q_ok, q_fail = gen_opq_ds(lg)

    lg.info("=======   STATISTICS   ========")
    if args.exps:
        lg.info(f"Stored experiment-lipid pairs: {e_ok} // failed: {e_fail}")
    if args.sims:
        lg.info(f"Stored simulation-lipid pairs: {s_ok} // failed: {s_fail}")
    if args.qual:
        lg.info(f"Stored paired simulation-lipid pairs: {q_ok} // failed: {q_fail}")

