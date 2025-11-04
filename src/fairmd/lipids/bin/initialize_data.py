#!/usr/bin/env python3
"""
Initializes FAIRMD Lipids data folder for further work and creates `databank_env.rc` file.


**Usage:**

.. code-block:: console

    fmdl_initialize_data [toy|stable|dev|-h]
    source databank_env.rc

- ``toy``: copies ToyData from the package data. This Databank contains only 4 very small systems
- ``stable``: download stable release from github.com/NMRlipids/BilayerData (even without git)
- ``dev``: clone NMRlipids/BilayerData (requires git)
"""

import importlib.util
import json
import os
import shutil
import subprocess
import sys
import zipfile
from io import BytesIO
from urllib.request import urlopen

DBREPO = "NMRlipids/BilayerData"
LOCDIRN = "BilayerData"
TOYDIRN = "ToyData"


def git_exists() -> bool:
    """Check if git is installed."""
    return shutil.which("git") is not None


def get_latest_release_tag() -> str:
    """Fetch the latest release tag from the GitHub API."""
    url = f"https://api.github.com/repos/{DBREPO}/releases/latest"
    with urlopen(url) as response:
        data = json.load(response)
        return data["tag_name"]


def clone_repo(tag: str) -> None:
    """Clone the repository at the specified tag."""
    if os.path.isdir(LOCDIRN):
        msg = f"Directory exists: {LOCDIRN}"
        raise FileExistsError(msg)
    subprocess.run(["git", "clone", "--branch", tag, f"https://github.com/{DBREPO}.git", LOCDIRN], check=True)


def download_and_extract_zip(tag: str) -> None:
    """Download and extract the ZIP file of the specified tag (release)."""
    unpack_name = DBREPO.split("/")[1] + "-" + tag
    for _name in [unpack_name, LOCDIRN]:
        if os.path.isdir(_name):
            msg = f"Directory exists: {_name}"
            raise FileExistsError(msg)
    url = f"https://github.com/{DBREPO}/archive/refs/tags/{tag}.zip"
    with urlopen(url) as response:
        with zipfile.ZipFile(BytesIO(response.read())) as z:
            z.extractall()
    os.rename(unpack_name, LOCDIRN)


def download_toy_folder(dest: str = TOYDIRN) -> None:
    """Download a folder from a package data."""
    pkgspec = importlib.util.find_spec("fairmd")
    src_path = os.path.join(pkgspec.submodule_search_locations[0], "lipids", "data", "ToyData")
    shutil.copytree(src_path, dest)


def initialize_data() -> None:
    """Code program logic."""
    print("== FAIRMD Lipids Data Initializer ==")
    # can be prg toy | prg stable | prg dev
    if len(sys.argv) != 2:
        sys.stderr.write("Usage: fmdl_initialize_data.py [toy|stable|dev]\n")
        sys.exit(1)

    prg = sys.argv[1].lower()
    if prg not in ["toy", "stable", "dev", "-h", "--help"]:
        sys.stderr.write("Invalid argument. Use 'toy', 'stable', or 'dev'.")
        sys.exit(1)
    sim_path = None
    data_path = None

    if prg == "toy":
        print("Downloading toy data ...")
        download_toy_folder()
        data_path = os.path.join(os.getcwd(), TOYDIRN)
        sim_path = os.path.join(data_path, "Simulations.1")
    elif prg == "stable":
        if git_exists():
            tag = get_latest_release_tag()
            print(f"Git found. Cloning {DBREPO} at {tag} ...")
            clone_repo(tag)
        else:
            tag = get_latest_release_tag()
            print(f"Git not found. Downloading ZIP for {DBREPO} at {tag} ...")
            download_and_extract_zip(tag)
            print("Done.")
        data_path = os.path.join(os.getcwd(), LOCDIRN)
    elif prg == "dev":
        if not git_exists():
            sys.stderr.write("Git is required for 'dev' option. Please install git.\n")
            sys.exit(1)
        print(f"Cloning {DBREPO} at main branch ...")
        clone_repo("main")
        data_path = os.path.join(os.getcwd(), LOCDIRN)
    elif prg in ["-h", "--help"]:
        print("""
Initializes FAIRMD Lipids data folder for further work and creates `databank_env.rc` file.

Usage:
    fmdl_initialize_data [toy|stable|dev|-h]
    source databank_env.rc

toy - copies ToyData from the package data. This Databank contains only 4 very small systems.
stable - download stable release from github.com/NMRlipids/BilayerData (even without git).
dev - clone NMRlipids/BilayerData (requires git).
""")

    # Write environment setup file
    with open("databank_env.rc", "w") as f:
        f.write(f"export FMDL_DATA_PATH={data_path}\n")
        if sim_path is not None:
            f.write(f"export FMDL_SIMU_PATH={sim_path}\n")

    print(f""""Data initialized into {data_path}. Please do

   $ source databank_env.rc

to set environment variables.""")


if __name__ == "__main__":
    initialize_data()
