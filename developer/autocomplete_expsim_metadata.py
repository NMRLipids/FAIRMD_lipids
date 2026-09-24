"""
Script for simulation and experiment Bioschemas autocomplete.

This script fills a ``bioschema_properties`` block in a simulation or experiment
``README.yaml`` so the record can be published as a `Bioschemas Dataset
<https://bioschemas.org/profiles/Dataset/1.0-RELEASE>`_. It is the dataset-level
counterpart of ``autocomplete_mol_metadata.py``, which does the same job for
molecules.

``name`` and ``description`` are composed from the record's own fields --
composition, temperature, hydration, ions, method -- and never taken from the
registry. A registry title names the paper or the Zenodo deposition, not the one
measurement or trajectory the record holds, so it describes the wrong thing and
does not tell sibling records apart: before this, 104 experiments shared 37
titles and 896 simulations shared 479. The fetched title is kept where it is
true, as ``isPartOf.name`` on the parent work, and fetched abstracts are dropped.

Every composed title ends in a bracketed tag that keeps near-identical records
apart -- the databank ``ID`` for a simulation, the first author and year for an
experiment. ``--check`` reports any two records that still share a title, and is
meant to be pointed at the whole databank, since only a full run can see a clash.

An experiment can carry two DOIs that mean different things, and the ``fields``
module decides once which is used for what: the ``ARTICLE_DOI`` is looked up,
because CrossRef holds the authors and the journal, and it becomes ``isPartOf``;
the ``DATA_DOI`` is what gets cited and what ``sameAs`` points at, because the
rule in ``docs/src/schemas/experiment_metadata.md`` is to cite the data. With
only one of them given, that one fills every role.

A plain ``DOI`` is the deprecated experiment spelling of those two, and is not
accepted by ``experiment_schema.json``. It is still read, so a record carrying
it keeps its citation, and both a run and ``--check`` name it wherever they find
one -- renaming it in the record is what clears the report. On a simulation the
same key is the Zenodo deposition and is not deprecated.

The remaining values are resolved from the record's DOI:

- DataCite -- Zenodo depositions, i.e. every simulation
- CrossRef -- journal articles, i.e. most experiments
- SPDX     -- the licence list, used as the licence controlled vocabulary

giving creators, dates, licence, publisher, citations and the parent work's
title. These are enriched from data already in the repository (composition,
force field, NMR/X-ray method, the analysis outputs present beside the README)
with terms from EDAM, CHMO and UO.

Properties that depend on where the databank is deployed -- ``identifier``,
``url``, ``@id``, ``@type``, ``@context``, ``dct:conformsTo`` and
``includedInDataCatalog`` -- are deliberately *not* written here; the web
frontend supplies them.

The script also normalises the surrounding record: dates are rewritten as quoted
ISO strings, and ``PUBLICATION`` is retired in favour of ``citation`` once its
content is safely represented there.

The work itself is not in this file. It lives in the ``expsim_metadata``
package beside it, one module per layer -- ``constants`` (with the ontology
terms and the databank's own names in ``constants.yaml``), ``helpers``,
``licenses``, ``registries``, ``fields``, ``descriptions``, ``bioschema``,
``records`` and ``checks``. This file is the command line over them.

.. note::
   This file is meant to be used by automated workflows.

   Unlike ``autocomplete_mol_metadata.py`` the file is **not** re-serialised. The
   generated block is appended, and only that block is rewritten in place, so
   hand-written comments elsewhere in the README survive -- most experiment
   files carry them. Runs are idempotent: a record whose content has not changed
   keeps its original ``_source.retrieved`` date.

   Upstream services answer with transient ``5xx`` errors from time to time, so
   requests are retried with exponential backoff that honours any
   ``Retry-After`` header. The retry budget can be overridden with the
   ``AUTOCOMPLETE_MAX_RETRIES`` environment variable (``0`` disables retries).
"""

import argparse
import os
import sys
from pathlib import Path

try:  # run as a script: developer/ is on sys.path
    from expsim_metadata.checks import check, duplicate_names
    from expsim_metadata.fields import data_root_of, molecule_names
    from expsim_metadata.licenses import load_spdx
    from expsim_metadata.records import process
except ImportError:  # imported as developer.autocomplete_expsim_metadata, e.g. by autodoc
    from developer.expsim_metadata.checks import check, duplicate_names
    from developer.expsim_metadata.fields import data_root_of, molecule_names
    from developer.expsim_metadata.licenses import load_spdx
    from developer.expsim_metadata.records import process


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("files", nargs="+", help="README.yaml paths to enrich")
    parser.add_argument("--cache", type=Path, default=None,
                        help="directory for cached API responses")
    parser.add_argument("--dry-run", action="store_true",
                        help="report what would change without writing")
    parser.add_argument("--check", action="store_true",
                        help="validate existing blocks instead of writing")
    parser.add_argument("--strict", action="store_true",
                        help="with --check, fail on a missing block rather than warn")
    args = parser.parse_args()

    paths = [Path(f) for f in args.files]
    missing = [p for p in paths if not p.is_file()]
    if missing:
        sys.exit(f"error: no such file: {missing[0]}")

    root = data_root_of(paths[0])
    cache_dir = args.cache or Path(os.environ.get("RUNNER_TEMP", root)) / ".cache" / "bioschema"
    spdx = load_spdx(cache_dir)
    if spdx is None:
        print("warning: SPDX licence list unavailable; licences will not resolve")

    if args.check:
        failed = set()
        for path in paths:
            for level, message in check(path, spdx, args.strict):
                print(message)
                if level == "ERROR":
                    failed.add(path)
        duplicates = duplicate_names(paths)
        for _, message in duplicates:
            print(message)
        clashes = [message for level, message in duplicates if level == "ERROR"]
        print(f"\n{len(paths) - len(failed)} of {len(paths)} records valid"
              + (f", {len(failed)} failing" if failed else "")
              + (f", {len(duplicates)} duplicated names" if duplicates else ""))
        sys.exit(1 if failed or clashes else 0)

    names = molecule_names(root)
    changed = sum(process(p, spdx, names, cache_dir, args.dry_run) for p in paths)
    print(f"\n{changed} of {len(paths)} records {'would be ' if args.dry_run else ''}updated")


if __name__ == "__main__":
    main()
