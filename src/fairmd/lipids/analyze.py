"""
fairmd.lipids.analyze contains top-level methods for properties' computing.

Module includes APL, TH, MAICoS, NMRPCA, and OP core analysis functions.
Code-intensive computations are performed inside corresponding modules.
Here, we gather only methods exported to the end-user if one wants to start
recomputing any of Databank's systems via python interface.
"""

import gc
import json
import os
import re
import socket
import subprocess
import urllib.request
from logging import Logger

import buildh
import MDAnalysis as mda
import numpy as np
from maicos.core.base import AnalysisCollection
from tqdm import tqdm

from fairmd.lipids import (
    FMDL_DATA_PATH,
    FMDL_SIMU_PATH,
    RCODE_COMPUTED,
    RCODE_ERROR,
    RCODE_SKIPPED,
)
from fairmd.lipids import analyze_nmrpca as nmrpca
from fairmd.lipids.core import System
from fairmd.lipids.databankio import download_resource_from_uri, resolve_download_file_url
from fairmd.lipids.databankLibrary import GetNlipids, getLipids, system2MDanalysisUniverse
from fairmd.lipids.databankop import find_OP
from fairmd.lipids.jsonEncoders import CompactJSONEncoder
from fairmd.lipids.maicos import (
    DensityPlanar,
    DielectricPlanar,
    DiporderPlanar,
    FormFactorPlanar,
    first_last_carbon,
    is_system_suitable_4_maicos,
    traj_centering_for_maicos,
)
from fairmd.lipids.settings import elements
from fairmd.lipids.settings.engines import get_struc_top_traj_fnames
from fairmd.lipids.settings.molecules import lipids_set


def computeNMRPCA(  # noqa: N802 (API)
    system: System,
    logger: Logger,
    *,
    recompute: bool = False,
) -> int:
    """Compute eq_times.json using NMR PCA analysis.

    :param system: System of the Databank
    :param logger: Logger object
    :param recompute: Delete previous apl.json and recompute it if True.
                      Defaults to False.

    :return: int success code (``RCODE_XXX``)
    """
    # getting data from databank and preprocessing them
    # Start Parser
    # TODO: 2test|    parser = Parser(FMDL_SIMU_PATH, readme, eq_time_fname, testTraj)
    try:
        parser = nmrpca.Parser(system, "eq_times.json")
        # Check trajectory
        print(parser._path)
        vpcode = parser.validate_path()
        print("ValidatePath code: ", vpcode)
    except Exception:
        logger.exception("Error initializing NMRPCA parser.")
        return RCODE_ERROR

    if vpcode > 0:
        logger.error("Some errors in analyze_nmrpca.py::Parser constructor.")
        return RCODE_ERROR

    if (
        (vpcode < 0 and not recompute)
        or ("WARNINGS" in system and system["WARNINGS"] is not None and "AMBIGUOUS_ATOMNAMES" in system["WARNINGS"])
        or ("WARNINGS" in system and system["WARNINGS"] is not None and "SKIP_EQTIMES" in system["WARNINGS"])
    ):
        return RCODE_SKIPPED

    # Check if TPR is defined and non-null or ""
    try:
        __ = system["TPR"][0][0]
        _ = __[0]
    except (KeyError, TypeError):
        logger.exception("TPR is required for NMRPCA analysis (%d})!", system["ID"])
        return RCODE_ERROR

    try:
        # Download files
        parser.download_traj()
        # Prepare trajectory
        parser.prepare_traj()
        # Concatenate trajectory
        parser.concatenate_traj()
        equilibration_times = {}
        # Iterate over trajectories for different lipids
        for lip, traj in parser.concatenated_trajs.items():
            print(f"Main: parsing lipid {lip}")
            # Creat PCA for trajectory
            pca_runner = nmrpca.PCA(traj[0], traj[1], traj[2], traj[3], parser.trj_len)
            # Run PCA
            data = pca_runner.PCA()
            print("Main: PCA: done")
            # Project trajectory on PC1
            pca_runner.get_proj(data)
            print("Main: Projections: done")
            # Calculate autocorrelations
            pca_runner.get_autocorrelations()
            print("Main: Autocorrelations: done")
            # Estimate equilibration time
            te2 = nmrpca.TimeEstimator(pca_runner.autocorrelation).calculate_time()
            equilibration_times[lip] = te2 / parser.trj_len
            print("Main: EQ time: done")

            print(te2 / parser.trj_len)

        parser.dump_data(equilibration_times)
        gc.collect()
    except Exception:
        logger.exception("Calculation of NMRPCA failed for %s", system["path"])
        return RCODE_ERROR
    else:
        return RCODE_COMPUTED


def computeAPL(  # noqa: N802 (API)
    system: System,
    logger: Logger,
    *,
    recompute: bool = False,
) -> int:
    """Generate apl.json analysis file for a system.

    Args:
        system (dict): one of systems of the Databank
        recompute (bool, optional): Delete previous apl.json and recompute it if True.
        Defaults to False.

    Returns
    -------
        int success code (``RCODE_XXX``)
    """
    # TODO: reading software and file path for simulations
    path = system["path"]

    # this is checking if area per lipid is already calculated for the systems
    outfilename = os.path.join(FMDL_SIMU_PATH, path, "apl.json")
    if os.path.isfile(outfilename):
        if recompute:
            os.unlink(outfilename)
        else:
            return RCODE_SKIPPED

    print("Analyzing: ", path)
    print("Will write into: ", outfilename)

    try:
        # calculates the total number of lipids
        n_lipid = GetNlipids(system)

        # makes MDAnalysis universe from the system. This also downloads the data if not
        # yet locally available
        u = system2MDanalysisUniverse(system)

        if u is None:
            print("Generation of MDAnalysis universe failed in folder", path)
            return RCODE_ERROR

        # this calculates the area per lipid as a function of time and stores it
        # in the databank
        apl = {}
        for _ts in tqdm(u.trajectory, desc="Scanning the trajectory"):
            if u.trajectory.time >= system["TIMELEFTOUT"] * 1000:
                dims = u.dimensions
                apl_frame = dims[0] * dims[1] * 2 / n_lipid
                apl[u.trajectory.time] = apl_frame

        with open(outfilename, "w") as f:
            json.dump(apl, f, cls=CompactJSONEncoder)
    except Exception:
        logger.exception("Calculation APL failed for %s ", system["path"])
        return RCODE_ERROR
    else:
        return RCODE_COMPUTED


def computeTH(  # noqa: N802 (API)
    system: System,
    logger: Logger,
    *,
    recompute: bool = False,
) -> int:
    thick_fn = os.path.join(FMDL_SIMU_PATH, system["path"], "thickness.json")
    if os.path.isfile(thick_fn) and not recompute:
        return RCODE_SKIPPED

    wat_dens_name = os.path.join(FMDL_SIMU_PATH, system["path"], "WaterDensity.json")
    lip_dens_name = os.path.join(FMDL_SIMU_PATH, system["path"], "LipidDensity.json")
    print(lip_dens_name)
    try:
        with open(wat_dens_name) as f:
            wat_dens = json.load(f)
        with open(lip_dens_name) as f:
            lip_dens = json.load(f)
        wd = np.array(wat_dens)
        ld = np.array(lip_dens)
        idx = np.argwhere(np.diff(np.sign(wd[:, 1] - ld[:, 1]))).flatten()
        if len(idx) < 2:
            print(
                "Dehydrated sample! Standard thickness rule doesn't work. Will extract boxZ.",
            )
            thickness = wd[-1, 0] - wd[0, 0]
        else:
            thickness = wd[idx[1], 0] - wd[idx[0], 0]
        with open(thick_fn, "w") as f:
            json.dump(thickness, f)
    except Exception:
        logger.exception("Calculation TH failed for %s.", system["path"])
        return RCODE_ERROR

    return RCODE_COMPUTED


# TODO: implement onlyLipid
def computeOP(  # noqa: N802 (API)
    system: System,
    logger: Logger,
    *,
    recompute: bool = False,
) -> int:
    """_summary_

    Args:
        system (dict): _description_
        recompute (bool, optional): _description_. Defaults to False.

    Returns
    -------
        int: _description_
    """
    path = system["path"]

    # Check if order parameters are calculated or something in the system prevents
    # order parameter calculation
    for key in system["COMPOSITION"]:
        outfilename = os.path.join(FMDL_SIMU_PATH, path, key + "OrderParameters.json")
        if os.path.isfile(outfilename) or (
            "WARNINGS" in system
            and type(system["WARNINGS"]) is dict
            and "AMBIGUOUS_ATOMNAMES" in system["WARNINGS"]
            and key in system["WARNINGS"]["AMBIGUOUS_ATOMNAMES"]
        ):
            file_found = True
        elif key in lipids_set:
            file_found = False
            break

    if file_found and not recompute:
        return RCODE_SKIPPED

    # If order parameters are not calculated and system ok, then continue to
    # calculating order parameters
    print("Analyzing: ", path)

    # Download the simulations and create MDAnalysis universe (the universe is not used
    # here but this script downloads the data)
    _ = system2MDanalysisUniverse(system)

    # Software and time for equilibration period
    software = system["SOFTWARE"]
    eq_time = float(system["TIMELEFTOUT"]) * 1000

    # Check if all or united atom simulation
    try:
        united_atom = system["UNITEDATOM_DICT"]
    except KeyError:
        united_atom = False

    # Check relevant warnings
    g3switch = (
        "WARNINGS" in system
        and type(system["WARNINGS"]) is dict
        and "GROMACS_VERSION" in system["WARNINGS"]
        and system["WARNINGS"]["GROMACS_VERSION"] == "gromacs3"
    )

    cur_path = os.path.join(FMDL_SIMU_PATH, path)

    try:
        struc_fname, top_fname, trj_fname = get_struc_top_traj_fnames(
            system,
            join_path=cur_path,
        )
    except (ValueError, KeyError):
        logger.exception("Error reading filenames from system dictionary.")
        return RCODE_ERROR

    try:
        # Set topology file names and make Gromacs trajectories whole
        if "gromacs" in software:
            if top_fname is None:
                raise ValueError("TPR is required for OP calculations!")
            xtcwhole = os.path.join(FMDL_SIMU_PATH, path, "whole.xtc")
            if not os.path.isfile(xtcwhole):
                try:
                    echo_proc = b"System\n"
                    if g3switch:
                        cmd_args = [
                            "trjconv",
                            "-f",
                            trj_fname,
                            "-s",
                            top_fname,
                            "-o",
                            xtcwhole,
                            "-pbc",
                            "mol",
                            "-b",
                            str(eq_time),
                        ]
                    else:
                        cmd_args = [
                            "gmx",
                            "trjconv",
                            "-f",
                            trj_fname,
                            "-s",
                            top_fname,
                            "-o",
                            xtcwhole,
                            "-pbc",
                            "mol",
                            "-b",
                            str(eq_time),
                        ]
                    if united_atom and system["TRAJECTORY_SIZE"] > 15e9:
                        cmd_args.extend(["-skip", "3"])
                    subprocess.run(cmd_args, input=echo_proc, check=True)
                except subprocess.CalledProcessError as e:
                    raise RuntimeError("trjconv exited with error (see above)") from e
        elif "openMM" in software or "NAMD" in software:
            if not os.path.isfile(struc_fname):
                pdb_url = resolve_download_file_url(system.get("DOI"), struc_fname)
                _ = urllib.request.urlretrieve(pdb_url, struc_fname)
        else:
            print(
                "Order parameter calculation for other than gromacs, openMM and NAMD are yet to be implemented.",
            )
            return RCODE_ERROR

        # Calculate order parameters
        # -----------------------------------
        # united-atom switch -> engine switch
        echo_proc = b"System\n"
        if united_atom:
            if "gromacs" not in software:
                raise ValueError("UNITED_ATOMS is supported only for GROMACS engine!")
            frame0struc = os.path.join(FMDL_SIMU_PATH, path, "frame0.gro")

            if g3switch:
                try:
                    subprocess.run(["editconf", "-f", top_fname, "-o", frame0struc], input=echo_proc, check=True)
                except subprocess.CalledProcessError as e:
                    raise RuntimeError("editconf exited with error (see above)") from e
            else:
                try:
                    if g3switch:
                        subprocess.run(
                            ["trjconv", "-f", xtcwhole, "-s", top_fname, "-dump", "0", "-o", frame0struc],
                            input=echo_proc,
                            check=True,
                        )
                    else:
                        subprocess.run(
                            ["gmx", "trjconv", "-f", xtcwhole, "-s", top_fname, "-dump", "0", "-o", frame0struc],
                            input=echo_proc,
                            check=True,
                        )
                except subprocess.CalledProcessError as e:
                    raise RuntimeError(
                        f"trjconv ({'trjconv' if g3switch else 'gmx trjconv'}) exited with error (see above)",
                    ) from e

            for key in system["UNITEDATOM_DICT"]:
                # construct order parameter definition file for CH bonds from
                # mapping file
                mapping_dict = system.content[key].mapping_dict

                def_fname = os.path.join(FMDL_SIMU_PATH, path, key + ".def")
                with open(def_fname, "w") as def_file:
                    previous_line = ""

                    regexp1_H = re.compile(r"M_[A-Z0-9]*C[0-9]*H[0-9]*_M")  # noqa: N806
                    regexp2_H = re.compile(r"M_G[0-9]*H[0-9]*_M")  # noqa: N806
                    regexp1_C = re.compile(r"M_[A-Z0-9]*C[0-9]*_M")  # noqa: N806
                    regexp2_C = re.compile(r"M_G[0-9]_M")  # noqa: N806

                    # Note that mapping_dict's keys must be ordered
                    # C11 H111 H112 C12 H121 H122 ...
                    # otherwize algorithm will fail
                    for mapping_key in mapping_dict:
                        if regexp1_C.search(mapping_key) or regexp2_C.search(mapping_key):
                            atom_c = [mapping_key, mapping_dict[mapping_key]["ATOMNAME"]]
                            atom_h = []
                        elif regexp1_H.search(mapping_key) or regexp2_H.search(mapping_key):
                            atom_h = [mapping_key, mapping_dict[mapping_key]["ATOMNAME"]]
                        else:
                            atom_c = []
                            atom_h = []

                        if atom_h:
                            assert atom_c
                            items = [atom_c[1], atom_h[1], atom_c[0], atom_h[0]]
                            def_line = items[2] + "&" + items[3] + " " + key + " " + items[0] + " " + items[1] + "\n"
                            if def_line != previous_line:
                                def_file.write(def_line)
                                previous_line = def_line

                # Add hydrogens to trajectory and calculate order parameters with buildH
                op_filepath = os.path.join(
                    FMDL_SIMU_PATH,
                    path,
                    key + "OrderParameters.dat",
                )

                lipid_json_file = [
                    os.path.join(
                        FMDL_DATA_PATH,
                        "lipid_json_buildH",
                        system["UNITEDATOM_DICT"][key] + ".json",
                    ),
                ]

                if not os.path.isfile(lipid_json_file[0]):
                    lipid_json_file = None

                print(system["UNITEDATOM_DICT"][key])
                buildh.launch(
                    coord_file=frame0struc,
                    def_file=def_fname,
                    lipid_type=system["UNITEDATOM_DICT"][key],
                    lipid_jsons=lipid_json_file,
                    traj_file=xtcwhole,
                    out_file=f"{op_filepath}.buildH",
                    ignore_CH3s=True,
                )

                with open(op_filepath + ".buildH") as op_file:
                    bh_lines = op_file.readlines()

                op_lines = []
                data = {}
                for line in bh_lines:
                    if "#" in line:
                        continue
                    line2 = (
                        line.split()[0].replace("&", " ")
                        + "  "
                        + line.split()[4]
                        + "  "
                        + line.split()[5]
                        + " "
                        + line.split()[6]
                        + "\n"
                    )
                    op_lines.append(line2)

                    op_name = line.split()[0].replace("&", " ")
                    op_values = [
                        float(line.split()[4]),
                        float(line.split()[5]),
                        float(line.split()[6]),
                    ]
                    data[str(op_name)] = []
                    data[str(op_name)].append(op_values)

                # write ascii file (TODO: DO WE NEED IT?)
                with open(op_filepath, "w") as outfile:
                    outfile.write("Atom     Average OP     OP stem\n")
                    outfile.writelines(op_lines)

                # write json
                outfile2 = os.path.join(
                    FMDL_SIMU_PATH,
                    path,
                    key + "OrderParameters.json",
                )
                with open(outfile2, "w") as f:
                    json.dump(data, f, cls=CompactJSONEncoder)

        # not united-atom cases
        else:
            if "gromacs" in software:
                gro = os.path.join(FMDL_SIMU_PATH, path, "conf.gro")

                print("\n Making gro file")
                if g3switch:
                    try:
                        subprocess.run(["editconf", "-f", top_fname, "-o", gro], input=echo_proc, check=True)
                    except subprocess.CalledProcessError as e:
                        raise RuntimeError("editconf exited with error (see above)") from e
                else:
                    try:
                        subprocess.run(
                            ["gmx", "trjconv", "-f", trj_fname, "-s", top_fname, "-dump", "0", "-o", gro],
                            input=echo_proc,
                            check=True,
                        )
                    except subprocess.CalledProcessError as e:
                        raise RuntimeError("trjconv exited with error (see above)") from e

            for key in system["COMPOSITION"]:
                if (
                    "WARNINGS" in system
                    and system["WARNINGS"] is not None
                    and "AMBIGUOUS_ATOMNAMES" in system["WARNINGS"]
                    and key in system["WARNINGS"]["AMBIGUOUS_ATOMNAMES"]
                ):
                    print(path, key)
                    print(
                        "Order parameters cannot be calculated if atom names are ambiguous.",
                    )
                    continue

                if key in lipids_set:
                    print("Calculating ", key, " order parameters")
                    resname = system["COMPOSITION"][key]["NAME"]
                    outfilename = os.path.join(
                        FMDL_SIMU_PATH,
                        path,
                        key + "OrderParameters.dat",
                    )
                    outfilename2 = os.path.join(
                        FMDL_SIMU_PATH,
                        path,
                        key + "OrderParameters.json",
                    )
                    if os.path.isfile(outfilename2):
                        print("Order parameter file already found")
                        continue
                    if "gromacs" in software:
                        try:
                            op_obj = find_OP(
                                system.content[key].mapping_dict,
                                top_fname,
                                xtcwhole,
                                resname,
                            )
                        except Exception as e:
                            logger.warning("We got this exception: %s", e)
                            logger.warning(
                                "But we will try rebuild the Universe from GROM if using tpr did not work!",
                            )
                            op_obj = find_OP(
                                system.content[key].mapping_dict,
                                gro,
                                xtcwhole,
                                resname,
                            )

                    if "openMM" in software or "NAMD" in software:
                        op_obj = find_OP(
                            system.content[key].mapping_dict,
                            struc_fname,
                            trj_fname,
                            resname,
                        )

                    data = {}

                    with open(outfilename, "w") as outfile:
                        outfile.write("Atom     Average OP     OP stem\n")

                        for _i, op in enumerate(op_obj):
                            (op.avg, op.std, op.stem) = op.get_avg_std_stem_OP
                            outfile.write(f"{op.name} {op.avg!s} {op.stem!s}\n")

                            data[str(op.name)] = []
                            data[str(op.name)].append(op.get_avg_std_stem_OP)

                    with open(outfilename2, "w") as f:
                        json.dump(data, f, cls=CompactJSONEncoder)

        print("Order parameters calculated and saved to ", path)

    except Exception:
        logger.exception("Calculation OP failed for %s.", system["path"])
        return RCODE_ERROR

    return RCODE_COMPUTED


def computeMAICOS(  # noqa: N802 (API)
    system: System,
    logger: Logger,
    *,
    ffonly: bool = True,
    recompute: bool = False,
) -> int:
    if not is_system_suitable_4_maicos(system):
        return RCODE_SKIPPED
    # otherwise continue
    software = system["SOFTWARE"]
    # download trajectory and gro files
    system_path = os.path.join(FMDL_SIMU_PATH, system["path"])
    doi = system.get("DOI")
    skip_downloading: bool = doi == "localhost"
    if skip_downloading:
        print("NOTE: The system with 'localhost' DOI should be downloaded by the user.")

    set_maicos_files = {
        "WaterDensity.json",
        "LipidDensity.json",
        "TotalDensity.json",
        "FormFactor.json",
    }

    if not ffonly:
        set_maicos_files |= {
            "LipidMassDensity.json",
            "WaterMassDensity.json",
            "TotalMassDensity.json",
            "DiporderWater.json",
            "Diporder2Water.json",
            "TotalChargeDensity.json",
            "WaterChargeDensity.json",
            "LipidChargeDensity.json",
            "DielectricLipid",  # Dielectric* files are treated
            "DielectricWater",  # slightly differently (see below)
            "DielectricTotal",
        }

    for file in set_maicos_files.copy():
        if "Dielectric" in file:
            if (
                os.path.isfile(os.path.join(system_path, file + "_par.json"))
                and os.path.isfile(os.path.join(system_path, file + "_perp.json"))
                and not recompute
            ):
                set_maicos_files.remove(file)
        elif os.path.isfile(os.path.join(system_path, file)) and not recompute:
            set_maicos_files.remove(file)

    try:
        struc, top, trj = get_struc_top_traj_fnames(system)
        trj_name = os.path.join(system_path, trj)
        struc_name = None if struc is None else os.path.join(system_path, struc)
        top_name = None if top is None else os.path.join(system_path, top)
    except Exception:
        logger.exception("Error getting structure/topology/trajectory filenames.")
        return RCODE_ERROR

    if top is None:
        # No topology => no charge inforamtion
        logger.info("Dielectric properties and charge densities are not accessible without topology.")
        for file in set_maicos_files.copy():
            if "Dielectric" in file or "Charge" in file:
                set_maicos_files.remove(file)
            if "Diporder" in file:
                # TODO: impute charges!
                set_maicos_files.remove(file)

    if not set_maicos_files:
        logger.info("All available MAICoS files are found. Skipping.")
        return RCODE_SKIPPED

    logger.info("Files to be computed: %s", "|".join(set_maicos_files))
    socket.setdefaulttimeout(15)

    def _not_skip_dwnld(f: str) -> bool:
        if skip_downloading and not os.path.isfile(f):
            msg = f"File [{trj_name}] should be downloaded by user"
            raise FileNotFoundError(msg)
        return not skip_downloading

    try:
        if _not_skip_dwnld(trj_name):
            trj_url = resolve_download_file_url(doi, trj)
            if not os.path.isfile(trj_name):
                print("Downloading trajectory with the size of ", system["TRAJECTORY_SIZE"], " to ", system["path"])
                _ = download_resource_from_uri(trj_url, trj_name)

        # make a function like this
        # TODO TPR should not be obligatory for GROMACS
        if "gromacs" in software:
            tpr_name = top_name
            if _not_skip_dwnld(tpr_name):
                tpr_url = resolve_download_file_url(doi, top)
                if not os.path.isfile(tpr_name):
                    _ = urllib.request.urlretrieve(tpr_url, tpr_name)
        elif "openMM" in software or "NAMD" in software:
            if _not_skip_dwnld(struc_name):
                pdb_url = resolve_download_file_url(doi, struc)
                if not os.path.isfile(struc_name):
                    _ = urllib.request.urlretrieve(pdb_url, struc_name)

        eq_time = float(system["TIMELEFTOUT"]) * 1000
        last_atom, g3_atom = first_last_carbon(system, logger)

        # Center around one lipid tail CH3 to guarantee all lipids in the same box
        if "gromacs" in system["SOFTWARE"]:
            # xtccentered
            xtccentered = traj_centering_for_maicos(
                system_path,
                trj_name,
                tpr_name,
                last_atom,
                g3_atom,
                eq_time,
                recompute=recompute,
            )
            u = mda.Universe(tpr_name, xtccentered)
        else:
            logger.warning("Centering for other than Gromacs is currently not implemented.")
            xtccentered = trj_name
            # it may not work w/o TPR if there are jumps over periodic boundary conditions in z-direction.
            u = mda.Universe(struc_name, xtccentered)

        # -- PHILIP code starts here --
        # We us a hardoced bin width
        bin_width = 0.3

        # introduce elements attribute (if it's empty)
        # and make initial guess (just in case)
        u.guess_TopologyAttrs(force_guess=["elements"])
        elements.guess_elements(system, u)

        # Adjust the group selection to be general for analysis
        water = u.select_atoms(f"resname {system['COMPOSITION']['SOL']['NAME']}")
        lipid = u.select_atoms(getLipids(system))
        # TODO: Maybe add group for ions and compute densities for them as well?

        # fixed `zmin` and `zmax` for profiles are deduced from smallest box dimension
        L_min = u.dimensions[2]  # noqa: N806 (PEP8)
        for ts in u.trajectory:
            L_min = min(L_min, ts.dimensions[2])

        # Skip unwrap/pack for speed - trajectories are already centered and whole
        base_options = {"unwrap": False, "bin_width": bin_width, "pack": False}
        zlim = {"zmin": -L_min / 2, "zmax": L_min / 2}
        dens_options = {**zlim, **base_options}

        spath = os.path.join(FMDL_SIMU_PATH, system["path"])

        request_analysis = {}

        if "TotalDensity.json" in set_maicos_files:
            # Density profiles
            dens_e_total = DensityPlanar(
                u.atoms,
                dens="electron",
                **dens_options,
                output=os.path.join(spath, "TotalDensity.json"),
            )
            request_analysis["TotalDensity.json"] = dens_e_total

        if "WaterDensity.json" in set_maicos_files:
            dens_e_water = DensityPlanar(
                water,
                dens="electron",
                **dens_options,
                output=os.path.join(spath, "WaterDensity.json"),
            )
            request_analysis["WaterDensity.json"] = dens_e_water

        if "LipidDensity.json" in set_maicos_files:
            dens_e_lipid = DensityPlanar(
                lipid,
                dens="electron",
                **dens_options,
                output=os.path.join(spath, "LipidDensity.json"),
            )
            request_analysis["LipidDensity.json"] = dens_e_lipid

        if "TotalMassDensity.json" in set_maicos_files:
            dens_m_total = DensityPlanar(
                u.atoms,
                dens="mass",
                **dens_options,
                output=os.path.join(spath, "TotalMassDensity.json"),
            )
            request_analysis["TotalMassDensity.json"] = dens_m_total

        if "WaterMassDensity.json" in set_maicos_files:
            dens_m_water = DensityPlanar(
                water,
                dens="mass",
                **dens_options,
                output=os.path.join(spath, "WaterMassDensity.json"),
            )
            request_analysis["WaterMassDensity.json"] = dens_m_water

        if "LipidMassDensity.json" in set_maicos_files:
            dens_m_lipid = DensityPlanar(
                lipid,
                dens="mass",
                **dens_options,
                output=os.path.join(spath, "LipidMassDensity.json"),
            )
            request_analysis["LipidMassDensity.json"] = dens_m_lipid

        if "TotalChargeDensity.json" in set_maicos_files:
            dens_c_total = DensityPlanar(
                u.atoms,
                dens="charge",
                **dens_options,
                output=os.path.join(spath, "TotalChargeDensity.json"),
            )
            request_analysis["TotalChargeDensity.json"] = dens_c_total

        if "WaterChargeDensity.json" in set_maicos_files:
            dens_c_water = DensityPlanar(
                water,
                dens="charge",
                **dens_options,
                output=os.path.join(spath, "WaterChargeDensity.json"),
            )
            request_analysis["WaterChargeDensity.json"] = dens_c_water

        if "LipidChargeDensity.json" in set_maicos_files:
            dens_c_lipid = DensityPlanar(
                lipid,
                dens="charge",
                **dens_options,
                output=os.path.join(spath, "LipidChargeDensity.json"),
            )
            request_analysis["LipidChargeDensity.json"] = dens_c_lipid

        # Form factor
        # Use `None` for `zmin`/`zmax` to respect changing cell size (vs fixed
        # zmin/zmax for density profiles)
        if "FormFactor.json" in set_maicos_files:
            form_factor = FormFactorPlanar(
                atomgroup=u.atoms,
                **base_options,
                zmin=None,
                zmax=None,
                output=os.path.join(spath, "FormFactor.json"),
            )
            request_analysis["FormFactor.json"] = form_factor

        # Water Orientation
        # TODO: Maybe also compute orientation for lipid heads/tails?
        if "DiporderWater.json" in set_maicos_files:
            cos_water = DiporderPlanar(
                water,
                order_parameter="cos_theta",
                **dens_options,
                output=os.path.join(spath, "DiporderWater.json"),
            )
            request_analysis["DiporderWater.json"] = cos_water

        if "Diporder2Water.json" in set_maicos_files:
            cos2_water = DiporderPlanar(
                water,
                order_parameter="cos_2_theta",
                **dens_options,
                output=os.path.join(spath, "Diporder2Water.json"),
            )
            request_analysis["Diporder2Water.json"] = cos2_water

        # Dielectric profiles
        if "DielectricTotal" in set_maicos_files:
            diel_total = DielectricPlanar(
                u.atoms,
                **base_options,
                output_prefix=os.path.join(spath, "DielectricTotal"),
            )
            request_analysis["DielectricTotal"] = diel_total

        if "DielectricWater" in set_maicos_files:
            diel_water = DielectricPlanar(
                water,
                **base_options,
                output_prefix=os.path.join(spath, "DielectricWater"),
            )
            request_analysis["DielectricWater"] = diel_water

        if "DielectricLipid" in set_maicos_files:
            diel_lipid = DielectricPlanar(
                lipid,
                **base_options,
                output_prefix=os.path.join(spath, "DielectricLipid"),
            )
            request_analysis["DielectricLipid"] = diel_lipid

        logger.info("We are ready to recompute %s", "|".join(request_analysis.keys()))
        if {"DielectricLipid", "DielectricWater", "DielectricTotal"}.intersection(request_analysis.keys()):
            logger.info("Check if dielectric profiles can be calculated (not possible for charged systems)")
            try:
                diel_total.run(stop=1)
            except (ValueError, UserWarning) as e:
                print(f"Dielectric profiles not available for this system: {e}")
                # create stub json-s to avoid recompute tries
                for dfile in ["DielectricTotal", "DielectricWater", "DielectricLipid"]:
                    with open(os.path.join(system_path, dfile + "_per.json"), "w") as f:
                        f.write("{}")
                    with open(os.path.join(system_path, dfile + "_par.json"), "w") as f:
                        f.write("{}")
                    logger.info(f"Created empty dielectric profile JSONs for {dfile}.")
                for k in ["DielectricTotal", "DielectricWater", "DielectricLipid"]:
                    request_analysis.pop(k, None)

        if request_analysis:
            collection = AnalysisCollection(*request_analysis.values())
            collection.run()

            for analysis in request_analysis.values():
                analysis.save()

    # finall catch
    except Exception:
        logger.exception("Calculation MAICOS failed for %s.", system["path"])
        return RCODE_ERROR
    else:
        return RCODE_COMPUTED
