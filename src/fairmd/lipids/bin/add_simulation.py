#!/usr/bin/env python3
r"""
Add new simulation into the FAIRMD Lipids.

The script adds a simulation into the Databank based on ``info.yaml`` file.

**Usage:**
    fmdl_add_simulation [-h] [-f FILE] [-d] [-n] [-w WORK_DIR] [--dry-run]

-h, --help             Show this help message and exit
-f FILE, --file=FILE   Input config file in yaml format.
-d, --debug            Enable debug logging output
-n, --no-cache         Always redownload repository files
-w WORK_DIR, --work-dir=WORK_DIR  Set custom temporary working directory \
        [not set = read from YAML]

**Returns** error codes:

- 0 - success
- 1 - input YAML parsing errors
- 2 - filesystem writting errors
- 3 - network accessing errors

"""

import argparse
import datetime
import importlib
import logging
import os
import pprint
import shutil
import subprocess
import sys
import tempfile
from copy import deepcopy
from urllib.error import HTTPError, URLError

import numpy as np
import pandas as pd
import yaml
from MDAnalysis import Universe

# import databank dictionaries
from fairmd.lipids import FMDL_SIMU_PATH
from fairmd.lipids.core import System, initialize_databank

# helpers
from fairmd.lipids.databankio import (
    calc_file_sha1_hash,
    create_databank_directories,
    download_resource_from_uri,
    resolve_download_file_url,
)
from fairmd.lipids.databankLibrary import lipids_set, molecules_set, parse_valid_config_settings
from fairmd.lipids.SchemaValidation.ValidateYAML import validate_info_dict
from fairmd.lipids.settings.engines import get_struc_top_traj_fnames, software_dict
from fairmd.lipids.settings.molecules import Lipid, NonLipid

pd.set_option("display.max_rows", 500)
pd.set_option("display.max_columns", 500)
pd.set_option("display.width", 1000)
pd.set_option("display.max_colwidth", 1000)


def add_simulation():
    info_template_path = os.path.join(
        importlib.util.find_spec("fairmd.lipids").submodule_search_locations[0],
        "SchemaValidation",
        "Schema",
        "info_template.yaml",
    )
    # parse input yaml file
    parser = argparse.ArgumentParser(
        prog="Add Simulation Script",
        description=f"""
Add new simulation into the FAIRMD Lipids. The script adds a
simulation into the Databank based on ``info.yaml`` file. You can
get template info-file from here:

{info_template_path}

Returns error codes:
- 0 - success
- 1 - input YAML parsing errors
- 2 - filesystem writting errors
- 3 - network accessing errors
""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-f",
        "--file",
        help="Input config file in yaml format.",
    )
    parser.add_argument(
        "-d",
        "--debug",
        help="enable debug logging output",
        action="store_true",
    )
    parser.add_argument(
        "-n",
        "--no-cache",
        help="always redownload repository files",
        action="store_true",
    )
    parser.add_argument(
        "-w",
        "--work-dir",
        help="set custom temporary working directory [not set = /tmp]",
        default=tempfile.gettempdir(),
    )
    parser.add_argument(
        "--dry-run",
        help="perform a dry-run download of the files with 50MB limit",
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

    input_path = os.path.normpath(args.file)

    # load input yaml file into empty dictionary
    info_yaml = {}

    # open input file for reading and writing
    with open(input_path) as yaml_file:
        info_yaml = yaml.load(
            yaml_file,
            Loader=yaml.FullLoader,  # noqa: S506
        )  # TODO may throw yaml.YAMLError

    errors = validate_info_dict(info_yaml)
    if errors:
        for error in errors:
            logger.exception(error)
        sys.exit(1)

    # Show the input read
    logger.debug(f"{os.linesep} Input read from {input_path} file:")
    pp = pprint.PrettyPrinter(width=41, compact=True)
    if logger.isEnabledFor(logging.DEBUG):
        pp.pprint(yaml.dump(info_yaml))

    # validate yaml entries and return updated sim dict
    try:
        sim_dict, files = parse_valid_config_settings(info_yaml)
    except KeyError:
        logger.exception("Missing entry key in yaml config, aborting..")
        sys.exit(1)
    except Exception as e:
        logger.exception(
            f"an '{type(e).__name__}' occured while performing validity check '{input_path}', script has been aborted",
        )
        sys.exit(1)
    else:
        logger.info(
            "all entries in simulation are understood and will be further processed",
        )
        logger.debug("valid sim entry keys:")
        pp = pprint.PrettyPrinter(width=41, compact=True)
        if logger.isEnabledFor(logging.DEBUG):
            pp.pprint(sim_dict)

    try:
        sim = System(sim_dict)
        # mapping files are registered here!
    except Exception as e:
        logger.exception(
            f"an '{type(e).__name__}' occured while processing dict->System '{input_path}', script has been aborted",
        )
        sys.exit(1)
    else:
        logger.info(f"System object is successfully created from '{input_path}' file")

    # Create temporary directory where to download files and analyze them
    dir_wrk = args.work_dir
    try:
        if args.no_cache:
            dir_tmp = tempfile.mkdtemp(prefix="tmp_6-", dir=dir_wrk)
        else:
            dir_tmp = os.path.join(dir_wrk, f"{sim['DOI'].split('/')[-1]}_download")
            os.makedirs(dir_tmp, exist_ok=True)
    except OSError:
        logger.exception("Couldn't create temporary working directory '%s'.", dir_tmp)
        sys.exit(2)
    logger.info(f"The data will be processed in directory path '{dir_tmp}'")

    # Check link status and download files
    try:
        download_links = []
        for fi in files:
            logger.info(f"Validating URL to file: {fi}..")
            _x = resolve_download_file_url(sim["DOI"], fi, validate_uri=True)
            download_links.append(_x)

        logger.info(f"Now downloading {len(files)} files ...")

        for url, fi in zip(download_links, files, strict=False):
            download_resource_from_uri(
                url,
                os.path.join(dir_tmp, fi),
                override_if_exists=args.no_cache,
                dry_run_mode=args.dry_run,
            )

        logger.info(f"Download of {len(files)} files was successful")

    except HTTPError as e:
        if e.code == 404:
            logger.exception(
                f"ressource not found on server '{e.url}' (404). Wrong DOI link or file name?",
            )
        else:
            logger.exception(f"HTTPError {e.code} while trying to download the file '{e.url}'")
        sys.exit(3)
    except URLError as e:
        logger.exception(
            f"couldn't resolve network adress: {e.reason}. Please check your internet connection.",
        )
        sys.exit(3)
    except Exception as e:
        logger.exception(
            f"'{type(e).__name__}' while attempting to download ressources, aborting",
        )
        sys.exit(3)

    # -- Calculate hash of downloaded files

    sim_hashes = deepcopy(sim)

    software_sim = software_dict[sim["SOFTWARE"].upper()]

    # list_containing the sha1 sums for all required files
    sha1_list_requied = []

    # Make empty dataframe with the desired columns
    df_files = pd.DataFrame(
        columns=["NAME", "TYPE", "REQUIRED", "SIZE_MB", "HASH"],
        dtype=object,
    )

    for key_sim, value_sim in sim_hashes.items():
        # double-checking keys
        try:
            entry_type = software_sim[key_sim]["TYPE"]
        except KeyError:
            if key_sim in ["SOFTWARE", "ID"]:
                continue
            else:
                # That shouldn't happen! Unexpected YAML-keys were checked by
                # parse_valid_config_settings before
                raise

        if "file" in entry_type:
            files_list = []
            is_required = software_dict[sim_hashes["SOFTWARE"].upper()][key_sim]["REQUIRED"]

            if not is_required and value_sim is None:
                continue  # skip not required NoneType (empty) file entries

            for file_provided in value_sim:
                file_name = os.path.join(dir_tmp, file_provided[0])
                logger.info(f"calculating sha1 hash of '{file_provided[0]}'...")
                file_hash = calc_file_sha1_hash(file_name)
                file_size_mb = f"{(os.path.getsize(file_name) / 1024 / 1024):.2f}"

                df_files = pd.concat(
                    [
                        df_files,
                        pd.DataFrame(
                            [
                                {
                                    "NAME": file_provided[0],
                                    "TYPE": key_sim,
                                    "REQUIRED": is_required,
                                    "SIZE_MB": file_size_mb,
                                    "HASH": file_hash,
                                },
                            ],
                        ),
                    ],
                    ignore_index=True,
                )
                files_list.append([file_provided[0], file_hash])

                # Find the keys of the required files to calculate the master_hash
                if is_required:
                    sha1_list_requied.append(file_hash)

                sim_hashes[key_sim] = files_list  # TODO Problematic

    logger.info(f"Summary of downloaded files: {os.linesep}")
    logger.info(df_files.to_string())

    # Calculate the hash of a file contaning the hashes of each of the required files
    # This should be always invariant as it will be used unique identifier for a
    # simualtion. Note order the hashes of the required files before calculating the
    # hash (That means that the required files cannot change)

    # Calculates numbers of lipid molecules in each leaflet. This is done by checking
    # on which side of the centre of mass the membrane each the centre of mass of a
    # lipid molecule is. If a lipid molecule is split so that headgroup and tails are
    # their own residues, the centre of mass of the headgroup is used in the
    # calculation.
    ####################################################################################

    logger.info(
        "Calculating the numbers of lipid molecules in each leaflet based on the "
        "center of mass of the membrane and lipids.",
    )
    logger.info(
        "If a lipid molecule is split to multiple residues, the centre of mass of the headgroup is used.",
    )

    try:
        struc, top, traj = get_struc_top_traj_fnames(sim, join_path=dir_tmp)
    except (ValueError, KeyError) as _:
        logger.exception("Some of fields required for Universe forming were not found.")
        sys.exit(1)
    except Exception as _:
        logger.exception("Unkonwn error during `get_3major_fnames..`")
        sys.exit(4)

    leaflet1 = 0  # total number of lipids in upper leaflet
    leaflet2 = 0  # total number of lipids in lower leaflet

    # try to generate zero-frame structure from top + trajectory
    fail_from_top = False
    gro = os.path.join(dir_tmp, "frame0.gro")  # structure regenerated from top
    try:
        logger.info(f"MDAnalysis tries to use {top} and {traj}")
        u = Universe(top, traj)
        u.atoms.write(gro, frames=u.trajectory[[0]])
    except Exception as e:
        logger.warning("%s: %s", e.__class__.__name__, e)
        fail_from_top = True

    # if previous fails then try the same from struc + trajectory
    if fail_from_top and struc is not None:
        try:
            logger.info(f"MDAnalysis tries to use {struc} and {traj}")
            u = Universe(struc, traj)
            u.atoms.write(gro, frames=u.trajectory[[0]])
        except Exception as _:
            logger.exception("Cannot initialize MDAnalysis using given structure file!")
            sys.exit(2)

    # if there is no struc and MDAnalysis doesn't start from TOP, then
    # GROMACS can try making struc from top!
    if fail_from_top and struc is None and sim["SOFTWARE"].upper() == "GROMACS":
        logger.info(
            "Now generating frame0.gro with Gromacs because MDAnalysis cannot read tpr version ...",
        )
        if "WARNINGS" in sim and sim["WARNINGS"] is not None and sim["WARNINGS"]["GROMACS_VERSION"] == "gromacs3":
            command = ["trjconv", "-s", top, "-f", traj, "-dump", "0", "-o", gro]
        else:
            command = ["gmx", "trjconv", "-s", top, "-f", traj, "-dump", "0", "-o", gro]
        logger.debug(f"executing 'echo System | {' '.join(command)}'")
        try:
            subprocess.run(command, input="System\n", text=True, check=True, capture_output=True)
        except subprocess.CalledProcessError as e:
            FAIL_MSG = f"Command 'echo System | {' '.join(command)}' failed with error: {e.stderr}"
            raise RuntimeError(FAIL_MSG) from e
        try:
            u = Universe(gro, traj)
            # write first frame into gro file
            u.atoms.write(gro, frames=u.trajectory[[0]])
        except Exception as e:
            logger.warning("%s: %s", e.__class__.__name__, e)
        struc = gro

    # if there is a topology and MDAnalysis reads it, we can use zero-frame
    # extraction as a structure!
    if not fail_from_top and struc is None:
        struc = gro

    # At last, we create universe from just a structure!
    logger.info(f"Making Universe from {struc}..")
    u0 = Universe(struc)

    lipids = []

    # select lipids
    for key_mol in lipids_set.names:
        logger.info(f"Calculating number of '{key_mol}' lipids")
        selection = ""
        if key_mol in sim["COMPOSITION"]:
            lip = Lipid(key_mol)
            m_file = sim["COMPOSITION"][key_mol]["MAPPING"]
            lip.register_mapping(m_file)
            for key in lip.mapping_dict:
                if "RESIDUE" in lip.mapping_dict[key]:
                    selection = (
                        selection
                        + "(resname "
                        + lip.mapping_dict[key]["RESIDUE"]
                        + " and name "
                        + lip.mapping_dict[key]["ATOMNAME"]
                        + ") or "
                    )
                else:
                    selection = "resname " + sim["COMPOSITION"][key_mol]["NAME"]
                    break
        selection = selection.removesuffix(" or ")
        molecules = u0.select_atoms(selection)
        if molecules.n_residues > 0:
            lipids.append(u0.select_atoms(selection))

    if not lipids:
        raise RuntimeError("No lipids were found in the composition!")
    # join all the selected the lipids together to make a selection of the entire
    # membrane and calculate the z component of the centre of mass of
    # the membrane
    membrane = u0.select_atoms("")
    R_membrane_z = 0
    if lipids:
        for i in range(len(lipids)):
            membrane = membrane + lipids[i]
        R_membrane_z = membrane.center_of_mass()[2]
    logger.info(f"Center of the mass of the membrane: {R_membrane_z!s}")

    # ---- number of each lipid per leaflet
    # TODO: remove code duplication block!
    for key_mol in lipids_set.names:
        leaflet1 = 0
        leaflet2 = 0

        selection: str = ""
        if key_mol in sim["COMPOSITION"]:
            lip = Lipid(key_mol)
            m_file = sim["COMPOSITION"][key_mol]["MAPPING"]
            lip.register_mapping(m_file)
            for key in lip.mapping_dict:
                if "RESIDUE" in lip.mapping_dict[key]:
                    selection = (
                        selection
                        + "resname "
                        + lip.mapping_dict[key]["RESIDUE"]
                        + " and name "
                        + lip.mapping_dict[key]["ATOMNAME"]
                        + " or "
                    )
                    break
                else:
                    selection = "resname " + sim["COMPOSITION"][key_mol]["NAME"]
                    break

        # if lipid was found then selection is not empty
        if selection != "":
            selection = selection.removesuffix(" or ")
            logger.debug(f"Selection: `{selection}`")
            molecules = u0.select_atoms(selection)
            logger.debug(
                "Resnames: %s | ResIDs: %s",
                ", ".join(molecules.residues.resnames),
                ", ".join(map(str, molecules.residues.resids)),
            )

            if molecules.n_residues > 0:
                for mol in molecules.residues:
                    R = mol.atoms.center_of_mass()

                    if R[2] - R_membrane_z > 0:
                        leaflet1 = leaflet1 + 1
                    elif R[2] - R_membrane_z < 0:
                        leaflet2 = leaflet2 + 1

        try:
            sim["COMPOSITION"][key_mol]["COUNT"] = [leaflet1, leaflet2]
        except KeyError:
            continue
        else:
            logger.info(f"Number of '{key_mol}' in upper leaflet: {leaflet1!s}")
            logger.info(f"Number of '{key_mol}' in lower leaflet: {leaflet2!s}")

    # ----- numbers of other molecules

    for key_mol in molecules_set.names:
        try:
            mol_name = sim["COMPOSITION"][key_mol]["NAME"]
        except KeyError:
            continue
        else:
            mol_number = u0.select_atoms("resname " + mol_name).n_residues
            sim["COMPOSITION"][key_mol]["COUNT"] = mol_number
            logger.info(
                f"Number of '{key_mol}': {sim['COMPOSITION'][key_mol]['COUNT']!s}",
            )

    # Anne: Read trajectory size and length

    sim["TRAJECTORY_SIZE"] = os.path.getsize(traj)

    n_frames = len(u.trajectory)
    timestep = u.trajectory.dt

    logger.info(f"Number of frames: {n_frames}")
    logger.info(f"Timestep: {timestep}")

    trj_length = n_frames * timestep

    sim["TRJLENGTH"] = trj_length

    # Read temperature from tpr
    if sim["SOFTWARE"].upper() == "GROMACS":
        file1 = os.path.join(dir_tmp, "tpr.txt")

        logger.info("Exporting information with gmx dump")
        if (
            "WARNINGS" in sim
            and sim["WARNINGS"] is not None
            and "GROMACS_VERSION" in sim["WARNINGS"]
            and sim["WARNINGS"]["GROMACS_VERSION"] == "gromacs3"
        ):
            command = ["gmxdump", "-s", top]
            TemperatureKey = "ref_t"
        else:
            command = ["gmx", "dump", "-s", top]
            TemperatureKey = "ref-t"
        try:
            result = subprocess.run(command, input="System\n", text=True, check=True, capture_output=True)
            with open(file1, "w") as f:
                f.write(result.stdout)
        except subprocess.CalledProcessError as e:
            FAIL_MSG = f"Command 'echo System | {' '.join(command)}' failed with error: {e.stderr}"
            raise RuntimeError(FAIL_MSG) from e

        with open(file1) as tpr_info:
            for line in tpr_info:
                if TemperatureKey in line:
                    sim["TEMPERATURE"] = float(line.split()[1])

    logger.info("Parameters read from input files:")
    logger.info(f"TEMPERATURE: {sim['TEMPERATURE']!s}")
    logger.info(f"LENGTH OF THE TRAJECTORY: {sim['TRJLENGTH']!s}")

    # Check that the number of atoms between data and README.yaml match

    natoms_trj = u.atoms.n_atoms

    number_of_atoms = 0
    for key_mol in sim["COMPOSITION"]:
        mol = Lipid(key_mol) if key_mol in lipids_set else NonLipid(key_mol)
        mol.register_mapping(sim["COMPOSITION"][key_mol]["MAPPING"])

        if sim.get("UNITEDATOM_DICT") and "SOL" not in key_mol:
            mapping_file_length = 0

            for key in mol.mapping_dict:
                if "H" in key:
                    continue
                else:
                    mapping_file_length += 1

        else:
            mapping_file_length = len(mol.mapping_dict)

        number_of_atoms += np.sum(sim["COMPOSITION"][key_mol]["COUNT"]) * mapping_file_length

    if number_of_atoms != natoms_trj:
        stop = input(
            f"Number of atoms in trajectory {natoms_trj} and README.yaml "
            f"{number_of_atoms} do no match. Check the mapping files and molecule"
            f" names. {os.linesep} If you know what you are doing, you can still "
            "continue the running the script. Do you want to (y/n)?",
        )
        if stop == "n":
            os._exit("Interrupted because atomnumbers did not match")
        if stop == "y":
            logger.warning(
                "Progressed even thought that atom numbers did not match. CHECK RESULTS MANUALLY!",
            )

    sim["NUMBER_OF_ATOMS"] = natoms_trj
    logger.info(f"Number of atoms in the system: {sim['NUMBER_OF_ATOMS']!s}")

    # ---- DATE OF RUNNING ----
    today = datetime.datetime.now().date().strftime("%d/%m/%Y")  # noqa: DTZ005
    sim["DATEOFRUNNING"] = today

    logger.info(f"Date of adding to the databank: {sim['DATEOFRUNNING']}")

    # Type of system is currently hard coded because only lipid bilayers are currently
    # added. When we go for other systems, this will be given by user.
    if "TYPEOFSYSTEM" not in list(sim.keys()):
        sim["TYPEOFSYSTEM"] = "lipid bilayer"

    # Determine inserting ID. We set -1 -2 -3 .. to new systems
    ss = initialize_databank()
    id_list = [s["ID"] for s in ss] + [0]
    sim["ID"] = min(id_list) - 1
    logger.info("Inserting ID: %s", str(sim["ID"]))
    del ss

    # dictionary saved in yaml format
    outfile_dict = os.path.join(dir_tmp, "README.yaml")
    with open(outfile_dict, "w") as f:
        yaml.dump(sim.readme, f, sort_keys=False, allow_unicode=True)
    logger.info(f"Databank README was saved to '{outfile_dict}'")

    # Try to create final directory
    try:
        directory_path = create_databank_directories(
            sim,
            sim_hashes,
            FMDL_SIMU_PATH,
            dry_run_mode=args.dry_run,
        )
    except NotImplementedError:
        logger.exception("[deprecated] Special error during directory creation (not implemented)")
        sys.exit(4)
    except OSError as e:
        logger.exception("couldn't create output directory: %s", e.args[1])
        sys.exit(2)

    logger.info("Databank entry will be registered into '%s'", directory_path)

    # copy previously downloaded files
    if not args.dry_run:
        logger.info("Copying files to the output directory [try hardlink for the traj.]...")
        try:
            os.link(
                traj,
                os.path.join(directory_path, os.path.basename(traj)),
            )
        except OSError:
            logger.warning(
                f"Could not hardlink trajectory file '{traj}' to the output directory. Copying instead.",
            )
            shutil.copyfile(traj, os.path.join(directory_path, os.path.basename(traj)))
        shutil.copyfile(
            top,
            os.path.join(directory_path, os.path.basename(top)),
        )
        shutil.copyfile(
            os.path.join(dir_tmp, "README.yaml"),
            os.path.join(directory_path, "README.yaml"),
        )

    logger.info("Script completed successfully!")


if __name__ == "__main__":
    add_simulation()
