#!/usr/bin/env python3

import argparse
import logging
import sys
from pathlib import Path

try:
    from . import _bootstrap
except ImportError:  # pragma: no cover - used when executed as a file
    import _bootstrap  # noqa: F401
from datasets.opqual import OPQDataStorer, gen_opq_ds
from datasets.pureop import ExpOPDataStorer, SimOPDataStorer, gen_op_from_exps, gen_op_from_sims


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser for dataset generation."""
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
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("."),
        help="Directory for generated HDF5 files (default: current directory)",
    )
    return parser


def configure_logger() -> logging.Logger:
    """Create the logger used by the publishing workflow."""
    lg = logging.getLogger("dataset-publishing")
    lg.setLevel(logging.INFO)
    if not lg.handlers:
        h_stderr = logging.StreamHandler(sys.stderr)
        h_stderr.setLevel(logging.INFO)
        lg.addHandler(h_stderr)
    return lg


def run(args: argparse.Namespace, log: logging.Logger | None = None) -> dict[str, tuple[int, int]]:
    """Generate the selected datasets and return per-dataset statistics."""
    if not args.exps and not args.sims and not args.qual:
        msg = "at least one of --exps, --sims, or --qual is required"
        raise ValueError(msg)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    log = log or configure_logger()
    results = {}
    if args.exps:
        results["exps"] = gen_op_from_exps(log, str(args.output_dir / ExpOPDataStorer.DEFAULT_EXP_HDFNAME))
    if args.sims:
        results["sims"] = gen_op_from_sims(log, str(args.output_dir / SimOPDataStorer.DEFAULT_SIMS_HDFNAME))
    if args.qual:
        results["qual"] = gen_opq_ds(log, str(args.output_dir / OPQDataStorer.DEFAULT_OPQ_HDFNAME))
    return results


def main(argv: list[str] | None = None) -> int:
    """Run the dataset-generation command-line workflow."""
    parser = build_parser()
    args = parser.parse_args(argv)

    if not args.exps and not args.sims and not args.qual:
        parser.print_usage()
        return 1

    results = run(args)
    lg = configure_logger()

    lg.info("=======   STATISTICS   ========")
    labels = {
        "exps": "experiment-lipid pairs",
        "sims": "simulation-lipid pairs",
        "qual": "paired simulation-lipid pairs",
    }
    for name, (ok, failed) in results.items():
        lg.info("Stored %s: %s // failed: %s", labels[name], ok, failed)
    return 0


if __name__ == "__main__":
    sys.exit(main())
