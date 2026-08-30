#!/usr/bin/env python3

import argparse
import logging
import sys

from .datasets.pureop import gen_op_from_exps, gen_op_from_sims

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="OP Dataset Generator",
        description="""
CI helper script for delivering dataset for Kaggle in HDF5 format.
Two DataFrames are stored for each Sim/Exp-lipid pair:
1. SMILES-aligned values for each atom of the molecule
2. Composition table of membrane part of the system""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--exps", action="store_true", help="Generate DS from experiments")
    parser.add_argument("--sims", action="store_true", help="Generate DS from simulations")
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

    lg.info("=======   STATISTICS   ========")
    if args.exps:
        lg.info(f"Stored experiment-lipid pairs: {e_ok} // failed: {e_fail}")
    if args.sims:
        lg.info(f"Stored simulation-lipid pairs: {s_ok} // failed: {s_fail}")

