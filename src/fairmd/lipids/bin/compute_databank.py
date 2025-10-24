#!/usr/bin/env python3
r"""
Compute properties for systems of the databank.

Script can compute selected set of properties for selected range of systems of the databank.
Properties available: APL, thickness, PCA relaxation times, MAICOS (density profiles,
water orientation, dielectric profiles, X-ray form factor).

**Usage:**

.. code-block:: console

    fmdl_compute_databank [--apl] [--nmrpca] [--maicos] [--thickness] \
        [--OP] [--range 0-1000] [--debug] [-h]


**Command line arguments:**

--range=ID_RANGE
    Range of system IDs to analyze, e.g. 0-1000. The interval `ID_RANGE` \
    can be opened at either end, e.g. ``*-50`` or ``25-*``. Use ``*-0`` to analyze all \
    just-added systems.
--apl        Compute APL (Area Per Lipid) for all systems.
--nmrpca     Compute NMR PCA for all systems.
--ff         Compute MAICOS electron density and form-factor for all systems.
--maicos     Compute MAICOS profiles for all systems.
--thickness  Compute Thickness for all systems.
--OP         Compute Order Parameter for all systems.
--debug      Enable debug logging.
--help       Show this help message and exit.

"""

import argparse
import logging
import sys

from fairmd.lipids.utils import run_analysis


def compute_databank():
    parser = argparse.ArgumentParser(
        prog="compute_databank.py Script",
        description="Compute computable properties",
    )
    parser.add_argument(
        "--apl",
        help="Compute APL (Area Per Lipid) for all systems",
        action="store_true",
    )
    parser.add_argument(
        "--nmrpca",
        help="Compute NMR PCA for all systems",
        action="store_true",
    )
    parser.add_argument(
        "--ff",
        help="Compute MAICOS electron density and form-factor for all systems",
        action="store_true",
    )
    parser.add_argument(
        "--maicos",
        help="Compute all MAICOS properties for all systems",
        action="store_true",
    )
    parser.add_argument(
        "--thickness",
        help="Compute Thickness for all systems",
        action="store_true",
    )
    parser.add_argument(
        "--OP",
        "--op",
        help="Compute Order Parameter for all systems",
        action="store_true",
    )
    parser.add_argument(
        "-r",
        "--range",
        help="Range of system IDs to analyze, e.g. 0-1000",
        type=str,
        default=None,
        metavar="ID_RANGE",
    )
    parser.add_argument(
        "-d",
        "--debug",
        help="Enable debug logging",
        action="store_true",
    )
    args = parser.parse_args()

    # configure logging
    logging_level = logging.DEBUG if args.debug else logging.INFO
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s]: %(message)s",
        datefmt="%I:%M:%S %p",
        level=logging_level,
    )
    logger = logging.getLogger()

    id_range = (None, None)
    if args.range:
        # range can be *-50 or 25-* or 25-38
        if args.range == "*":
            id_range = (None, None)
        elif args.range.startswith("*-"):
            try:
                id_range = (None, int(args.range[2:]))
            except ValueError:
                logger.exception("Invalid ID range format. Use 'start-end' format.")
                sys.exit(1)
        elif args.range.endswith("-*"):
            try:
                id_range = (int(args.range[:-2]), None)
            except ValueError:
                logger.exception("Invalid ID range format. Use 'start-end' format.")
                sys.exit(1)
        else:
            # assume format is start-end
            try:
                id_range = tuple(map(int, args.range.split("-")))
            except ValueError:
                logger.exception("Invalid ID range format. Use 'start-end' format.")
                sys.exit(1)

    if args.apl:
        logger.info("Computing APL (Area Per Lipid) for all systems")
        from fairmd.lipids.analyze import computeAPL

        run_analysis(computeAPL, logger, id_range=id_range)
    if args.nmrpca:
        logger.info("Computing NMR PCA for all systems")
        from fairmd.lipids.analyze import computeNMRPCA

        run_analysis(computeNMRPCA, logger, id_range=id_range)

    if args.ff and args.maicos:
        pass

    if args.ff and not args.maicos:
        logger.info("Computing MAICoS electron density and form-factor for all systems")
        from fairmd.lipids.analyze import computeMAICOS

        run_analysis(computeMAICOS, logger, id_range=id_range)

    if args.maicos:
        logger.info("Computing MAICOS for all systems")
        from fairmd.lipids.analyze import computeMAICOS

        def compute_all_maicos_props(s, ell):
            return computeMAICOS(s, ell, ffonly=False)

        run_analysis(compute_all_maicos_props, logger, id_range=id_range)

    if args.thickness:
        logger.info("Computing Thickness for all systems")
        from fairmd.lipids.analyze import computeTH

        run_analysis(computeTH, logger, id_range=id_range)
    if args.OP:
        logger.info("Computing Order Parameter for all systems")
        from fairmd.lipids.analyze import computeOP

        run_analysis(computeOP, logger, id_range=id_range)


if __name__ == "__main__":
    compute_databank()
