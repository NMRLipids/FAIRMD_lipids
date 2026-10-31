#!/usr/bin/env python3

"""
Script for processing lipid databank systems and saving analysis results to HDF5.

This script filters systems by water-to-lipid ratio, calculates density
profiles using maicos, and exports the data for machine learning training.
"""

import argparse
import logging
import os
import re

import h5py
import numpy as np

from fairmd.lipids import FMDL_MAICOS_NCORES, FMDL_SIMU_PATH
from fairmd.lipids.analib.maicos import (
    DensityPlanar,
    FormFactorPlanar,
    first_last_carbon,
    is_system_suitable_4_maicos,
    traj_centering_for_maicos_gromacs,
    traj_centering_for_maicos_mda,
    traj_centering_for_maicos_mda_parallel,
)
from fairmd.lipids.api import UniverseConstructor, get_mean_ApL, get_thickness
from fairmd.lipids.auxiliary import mollib
from fairmd.lipids.core import System, initialize_databank


def is_suitable(system: System) -> bool:
    """
    Check if a simulation system is suitable for analysis with maicos.

    Checks for the presence of a 'WARNINGS' dictionary in the system metadata
    and verifies compatibility using the internal fairmd suitability check.

    Args:
        system (dict): The system dictionary from the databank.

    Returns:
        bool: True if the system is suitable for analysis, False otherwise.
    """
    flag = True
    if "WARNINGS" in system and isinstance(system["WARNINGS"], dict):
        flag = False
    if not is_system_suitable_4_maicos(system):
        print(f"system {system} not suitable for maicos")
        flag = False
    return flag


def get_scalar_properties(system: System) -> tuple[float, float, bool]:
    """
    Retrieve scalar physical properties (ApL and thickness) for a given system.

    Attempts to load the mean Area per Lipid (ApL) and bilayer thickness.
    If a property cannot be loaded, it assigns a default value of -1 and
    logs a warning.

    Args:
        system (dict): The system dictionary from the databank.

    Returns:
        tuple: A tuple containing:
            - ApL (float): Average area per lipid (Å²).
            - thickness (float): Bilayer thickness (Å).
            - no_error_flag (bool): True if all properties were loaded successfully.
    """
    no_error_flag = True
    try:
        apl = get_mean_ApL(system)
    except Exception as e:
        print(f"System {system} - can't load ApL: {e}")
        no_error_flag = False
        apl = -1

    try:
        thickness = get_thickness(system)
    except Exception as e:
        print(f"System {system} - can't load thickness: {e}")
        no_error_flag = False
        thickness = -1
    return apl, thickness, no_error_flag


def compute_pp_thickness(u, logger=None) -> float:
    """
    Compute phosphorus-phosphorus bilayer thickness, averaged over the trajectory.

    For each frame, phosphorus atoms are split into the two leaflets by their
    z position relative to the box center (the bilayer midplane established by
    `center_trajectory`), each leaflet's average z position is taken, and the
    distance between the two leaflet averages gives the frame's P-P thickness.

    Args:
        u (MDAnalysis.Universe): Centered universe with elements guessed
            (so phosphorus atoms can be selected by element).
        logger (logging.Logger, optional): Logger for warnings.

    Returns:
        float: Mean P-P thickness (Å) over the trajectory, or NaN if no
            phosphorus atoms are present or leaflets cannot be resolved.
    """
    p_atoms = u.select_atoms("element P")
    if len(p_atoms) == 0:
        if logger:
            logger.warning("No phosphorus atoms found, skipping P-P thickness")
        return np.nan

    frame_thickness = []
    for _ts in u.trajectory:
        box_center_z = u.dimensions[2] / 2
        z = p_atoms.positions[:, 2]
        upper = z[z >= box_center_z]
        lower = z[z < box_center_z]
        if len(upper) == 0 or len(lower) == 0:
            continue
        frame_thickness.append(abs(upper.mean() - lower.mean()))

    if not frame_thickness:
        if logger:
            logger.warning("Could not split phosphorus atoms into two leaflets")
        return np.nan
    return float(np.mean(frame_thickness))


def compute_peak_to_peak_thickness(bin_pos, profile) -> float:
    """
    Distance between the two leaflet peaks of a time-averaged planar density profile.

    Splits the profile at the bilayer midplane (``bin_pos == 0``, matching MAICoS'
    planar-analysis convention) and finds the position of the maximum value on
    each side, then returns the distance between the two peak positions.

    Args:
        bin_pos (array-like): Bin positions along the membrane normal (Å),
            centered on the bilayer midplane.
        profile (array-like): Density profile values matching ``bin_pos``.

    Returns:
        float: Distance between the lower- and upper-leaflet peaks (Å), or NaN
            if one side of the profile has no bins.
    """
    bin_pos = np.asarray(bin_pos)
    profile = np.asarray(profile)

    lower_mask = bin_pos < 0
    upper_mask = bin_pos >= 0
    if not lower_mask.any() or not upper_mask.any():
        return np.nan

    lower_peak_pos = bin_pos[lower_mask][np.argmax(profile[lower_mask])]
    upper_peak_pos = bin_pos[upper_mask][np.argmax(profile[upper_mask])]
    return float(upper_peak_pos - lower_peak_pos)


def center_trajectory(
    system: System,
    uc: UniverseConstructor,
    last_atom,
    g3_atom,
    eq_time,
    logger,
    *,
    recompute: bool = False,
):
    """
    Centers the simulation trajectory for analysis, handling different software backends.

    Coordinates trajectory centering using either Gromacs commands or MDAnalysis
    (sequential or parallel) based on the simulation metadata and environment
    configuration, then loads the centered trajectory into a Universe with
    elements guessed.

    Args:
        system (System): The system whose trajectory is being centered.
        uc (UniverseConstructor): Object containing simulation path and topology info.
        last_atom (int/str): Index or name of the last atom for centering reference.
        g3_atom (int/str): Index or name of the glycerol 3 atom for orientation.
        eq_time (float): Equilibration time to skip in milliseconds.
        logger (logging.Logger): Logger instance for status and error reporting.
        recompute (bool): If True, force recomputation even if a centered
            trajectory file already exists.

    Returns:
        MDAnalysis.Universe: The centered universe, loaded with the new
            trajectory and elements guessed.
    """
    u = uc.build_universe()
    spath = os.path.join(FMDL_SIMU_PATH, system["path"])
    if "gromacs" in system["SOFTWARE"]:
        traj_centered = traj_centering_for_maicos_gromacs(
            spath,
            tpr_name=uc.paths["top"],
            trj_name=uc.paths["traj"],
            last_atom=last_atom,
            g3_atom=g3_atom,
            eq_time=eq_time,
            recompute=recompute,
        )
    elif FMDL_MAICOS_NCORES != 1:
        try:
            n_jobs = FMDL_MAICOS_NCORES if FMDL_MAICOS_NCORES is not None else -1
            logger.info(f"Using parallel trajectory centering (n_jobs={n_jobs})")
            traj_centered = traj_centering_for_maicos_mda_parallel(
                u,
                spath,
                last_atom,
                eq_time,
                n_jobs=n_jobs,
                recompute=recompute,
                logger=logger,
                show_progress=True,
            )
        except ImportError:
            logger.warning("joblib not available, falling back to sequential centering")
            traj_centered = traj_centering_for_maicos_mda(
                u,
                spath,
                last_atom,
                eq_time,
                recompute=recompute,
                logger=logger,
            )
    else:
        logger.info("Using sequential trajectory centering (FMDL_MAICOS_NCORES=1)")
        traj_centered = traj_centering_for_maicos_mda(
            u,
            spath,
            last_atom,
            eq_time,
            recompute=recompute,
            logger=logger,
        )
    u.load_new(traj_centered, format="XTC")
    u.guess_TopologyAttrs(force_guess=["elements"])
    mollib.guess_elements(system, u)
    return u


def separate_lipid_atoms(mapping_dict):
    """
    Group lipid atom names into headgroup, tail, and backbone fragments.

    Parses a mapping dictionary (usually from a Lipid class instance) to
    categorize atoms based on their structural fragment.

    Args:
        mapping_dict (dict): Dictionary mapping atom IDs to names and fragments.

    Fragment labels for lipids with more than one copy of a fragment (e.g. cardiolipin's
    two glycerol backbones/tail pairs) carry a trailing index, e.g. "sn-1 2" or
    "glycerol backbone 1". The trailing " <number>" is stripped before matching, so
    those atoms are grouped with their un-indexed counterparts.

    Returns:
        tuple: A triplet of space-separated strings (head_atoms, tail_atoms, backbone_atoms)
            containing the atom names for each respective fragment.

    Raises:
        ValueError: If an atom's fragment label (after stripping any trailing index)
            doesn't match a known fragment category.
    """
    head_atoms = ""
    tail_atoms = ""
    backbone_atoms = ""
    for atom in mapping_dict:
        fragment = mapping_dict[atom]["FRAGMENT"]
        atom_name = mapping_dict[atom]["ATOMNAME"]
        base_fragment = re.sub(r"\s+\d+$", "", fragment)
        if base_fragment in ("headgroup", "glucose"):
            head_atoms += atom_name + " "
        elif base_fragment in ("sn-1", "sn-2", "sn-4", "tail"):
            tail_atoms += atom_name + " "
        elif base_fragment in ("glycerol backbone", "backbone"):
            backbone_atoms += atom_name + " "
        else:
            msg = f"Invalid atom - {atom} - {atom_name} - {fragment}"
            raise ValueError(msg)
    return (head_atoms, tail_atoms, backbone_atoms)


def create_fragment_selectors(system: System):
    """
    Create MDAnalysis atom selection strings for lipid fragments across a system.

    Iterates through a system's lipid objects (already registered against the
    mapping file this specific system actually uses, see ``System._initialize_content``)
    and constructs, per fragment, a selection combining every lipid's fragment
    atoms - each scoped to that lipid's own resname, joined with 'or' - so the
    resulting group is all lipids' headgroups (etc.) together, without atom-name
    collisions across lipid species pulling in the wrong atoms.

    Args:
        system (System): The system whose lipids' fragments are being selected.

    Returns:
        list of str: A list containing three selection strings in the order:
            [head_selector, tail_selector, backbone_selector].
    """
    head_clauses, tail_clauses, backbone_clauses = [], [], []
    for lipid_key, lipid_class in system.lipids.items():
        resname = system["COMPOSITION"][lipid_key]["NAME"]
        mapping_dict = lipid_class.mapping_dict
        head_atoms, tail_atoms, backbone_atoms = separate_lipid_atoms(mapping_dict)

        if head_atoms.strip():
            head_clauses.append(f"(name {head_atoms}and resname {resname})")
        if tail_atoms.strip():
            tail_clauses.append(f"(name {tail_atoms}and resname {resname})")
        if backbone_atoms.strip():
            backbone_clauses.append(f"(name {backbone_atoms}and resname {resname})")

    return [" or ".join(head_clauses), " or ".join(tail_clauses), " or ".join(backbone_clauses)]


def lipids_missing_fragments(lipids) -> bool:
    """
    Check whether any lipid in lipids lacks head, tail, or backbone atoms.

    Args:
        lipids (Iterable[Lipid]): Lipid molecule objects present in the system,
            e.g. ``system.lipids.values()``, already registered against the
            mapping file this specific system actually uses.

    Returns:
        bool: True if any lipid has an empty head, tail, or backbone selection.
    """
    for lipid_class in lipids:
        head_atoms, tail_atoms, backbone_atoms = separate_lipid_atoms(lipid_class.mapping_dict)
        if not head_atoms.strip() or not tail_atoms.strip() or not backbone_atoms.strip():
            return True
    return False


class HDF5LipidWriter:
    """
    Save lipid simulation analysis results into HDF5 format.

    This class manages the hierarchical storage of form factors, density
    profiles, and scalar properties. It is designed for long-running
    processes by opening and flushing to the file for every system processed.
    """

    def __init__(self, filename: str, *, overwrite_file: bool = False) -> None:
        """
        Initialize the writer and optionally clears the existing file.

        Args:
            filename (str): Path to the output .h5 file.
            overwrite_file (bool): If True, deletes the existing file on initialization.
        """
        self.filename = filename

        if overwrite_file and os.path.exists(self.filename):
            print(f"Clearing existing file: {self.filename}")
            os.remove(self.filename)

    def has_system(self, sys_id: str) -> bool:
        """
        Check whether a system's results are already stored in the HDF5 file.

        Args:
            sys_id (str): System ID as used for the top-level HDF5 group key.

        Returns:
            bool: True if the file exists and already contains this system.
        """
        if not os.path.exists(self.filename):
            return False
        with h5py.File(self.filename, "r") as f:
            return sys_id in f

    def save_system(self, system, scalar_data, form_factor, total_dens, mol_densities, frag_densities):
        """
        Save a single system's results to the HDF5 file.

        Organizes data into groups for metadata, axes, form factors, and
        various electron density profiles (total, fragment-based, and molecule-based).
        If the system ID already exists, it is overwritten to ensure data integrity.

        Args:
            system (dict): System metadata from the databank.
            scalar_data (dict): Scalar system properties, stored verbatim as group
                attributes (name -> value). Expected to contain 'ApL', 'thickness'
                (water/lipid density intersection), 'thickness_pp'
                (phosphorus-phosphorus), 'thickness_headgroup_peaks', and
                'thickness_totaldensity_peaks'.
            form_factor (tuple): (q_pos, profile, dprofile) for the form factor.
            total_dens (tuple): (r_pos, profile, dprofile) for total electron density.
            mol_densities (list of tuples): List of (name, profile, dprofile) for each
                molecule type, where name is the resname used to select it.
            frag_densities (list of tuples): List of (profile, dprofile) for lipid fragments.
        """
        with h5py.File(self.filename, "a") as f:
            sys_id = str(system["ID"])

            if sys_id in f:
                print(f"Warning: System {sys_id} already exists in HDF5. Overwriting.")
                del f[sys_id]

            grp = f.create_group(sys_id)

            grp.attrs["path"] = system.get("path", "")
            for key, value in scalar_data.items():
                grp.attrs[key] = value if value is not None else np.nan

            axis_grp = grp.create_group("axis")
            self._write_dataset(axis_grp, "q_pos", form_factor[0])
            self._write_dataset(axis_grp, "r_pos", total_dens[0])

            ff_grp = grp.create_group("form_factor")
            self._write_dataset(ff_grp, "profile", form_factor[1])
            self._write_dataset(ff_grp, "dprofile", form_factor[2])

            td_grp = grp.create_group("density_total")
            self._write_dataset(td_grp, "profile", total_dens[1])
            self._write_dataset(td_grp, "dprofile", total_dens[2])

            frag_labels = ["head", "tail", "backbone"]
            frag_grp = grp.create_group("density_fragments")
            for i, (profile, dprofile) in enumerate(frag_densities):
                label = frag_labels[i] if i < len(frag_labels) else f"frag_{i}"
                sub_grp = frag_grp.create_group(label)
                self._write_dataset(sub_grp, "profile", profile)
                self._write_dataset(sub_grp, "dprofile", dprofile)

            mol_grp = grp.create_group("density_molecules")
            for i, (name, profile, dprofile) in enumerate(mol_densities):
                sub_grp = mol_grp.create_group(f"mol_{i}")
                sub_grp.attrs["name"] = name
                self._write_dataset(sub_grp, "profile", profile)
                self._write_dataset(sub_grp, "dprofile", dprofile)

    def _write_dataset(self, group, name, data):
        """
        Write a NumPy array to an HDF5 group with compression.

        Args:
            group (h5py.Group): The parent group to write into.
            name (str): The name of the dataset.
            data (array-like): The numerical data to save.
        """
        if data is not None:
            group.create_dataset(name, data=np.array(data), compression="gzip", compression_opts=4)


def recompute_extended_ff_dataset(
    h5fpath: str, *,
    hydration_threshold: int = 20,
    recompute_centering: bool = True,
    small_trajs_only: bool = False,
    overwrite: bool = False,
    force_ids: set[str] | None = None,
) -> None:
    """
    Recompute the extended form factor dataset for all systems in the databank.

    This function iterates through all systems, checks their suitability, and
    processes them to extract form factors and density profiles. The results
    are saved into an HDF5 file using the HDF5LipidWriter class. Systems that
    do not meet the criteria (e.g., water-to-lipid ratio) are skipped. Systems
    already present in the output file are skipped as well, unless
    ``overwrite`` is set, or their ID is listed in ``force_ids``.
    """
    systems = initialize_databank()
    logger = logging.getLogger(__name__)
    writer = HDF5LipidWriter(h5fpath)
    force_ids = force_ids or set()

    count = 0
    print(f"Number of systems: {len(systems)}")
    for system in systems:
        print(f"Calculating system ID {system['ID']}, hash {system['path']}")

        if not overwrite and str(system["ID"]) not in force_ids and writer.has_system(str(system["ID"])):
            continue

        if system["TRAJECTORY_SIZE"] > 10**8 and small_trajs_only:  # For testing purpouses
            continue

        if not is_suitable(system):
            continue

        ApL, thickness, flag = get_scalar_properties(system)
        scalar_info = {"ApL": ApL, "thickness": thickness}

        if system.get_hydration(basis="number") < hydration_threshold:
            continue

        try:
            uc = UniverseConstructor(system)
            uc.download_mddata()
        except Exception as e:
            print(f"System {system} - can't build/download universe: {e}")
            continue

        eq_time = float(system["TIMELEFTOUT"]) * 1000
        last_atom, g3_atom = first_last_carbon(system, logger)

        u = center_trajectory(system, uc, last_atom, g3_atom, eq_time, logger, recompute=recompute_centering)

        print("Calculating P-P thickness")
        scalar_info["thickness_pp"] = compute_pp_thickness(u, logger=logger)

        bin_width = 0.3

        L_min = u.dimensions[2]
        for ts in u.trajectory:
            L_min = min(L_min, ts.dimensions[2])

        base_options = {"unwrap": False, "bin_width": bin_width, "pack": False}
        zlim = {"zmin": -L_min / 2, "zmax": L_min / 2}
        dens_options = {**zlim, **base_options}

        print("Calculating form factor")
        form_factor = FormFactorPlanar(
            atomgroup=u.atoms,
            **base_options,
            zmin=None,
            zmax=None,
        ).run()
        ff = (form_factor.results.bin_pos, form_factor.results.profile, form_factor.results.dprofile)

        print("Calculating total density")
        dens_total_runner = DensityPlanar(
            u.atoms,
            dens="electron",
            **dens_options,
        ).run()
        dens_total = (
            dens_total_runner.results.bin_pos,
            dens_total_runner.results.profile,
            dens_total_runner.results.dprofile,
        )

        molecule_names = [system["COMPOSITION"][molkey]["NAME"] for molkey in system.content]
        dens_molecule = []
        for name in molecule_names:
            print(f"Calculating {name} density")
            molecule_group = u.select_atoms(f"resname {name}")
            dens_molecule_runner = DensityPlanar(
                molecule_group,
                dens="electron",
                **dens_options,
            ).run()
            dens = (name, dens_molecule_runner.results.profile, dens_molecule_runner.results.dprofile)
            dens_molecule.append(dens)

        try:
            fragment_selectors = create_fragment_selectors(system)
        except ValueError as e:
            print(f"System {system} - skipping, {e}")
            continue
        dens_fragment = []
        frag_labels = ["head", "tail", "backbone"]
        for i, selector in enumerate(fragment_selectors):
            print(f"Calculating {frag_labels[i]} density")
            print(selector)
            fragment_group = u.select_atoms(selector)
            dens_fragment_runner = DensityPlanar(
                fragment_group,
                dens="electron",
                **dens_options,
            ).run()
            dens = (dens_fragment_runner.results.profile, dens_fragment_runner.results.dprofile)
            dens_fragment.append(dens)

        print("Calculating headgroup peak-to-peak thickness")
        head_profile = dens_fragment[frag_labels.index("head")][0]
        scalar_info["thickness_headgroup_peaks"] = compute_peak_to_peak_thickness(dens_total[0], head_profile)

        print("Calculating total density peak-to-peak thickness")
        scalar_info["thickness_totaldensity_peaks"] = compute_peak_to_peak_thickness(dens_total[0], dens_total[1])

        writer.save_system(
            system=system,
            scalar_data=scalar_info,
            form_factor=ff,
            total_dens=dens_total,
            mol_densities=dens_molecule,
            frag_densities=dens_fragment,
        )

        count += 1

    print(f"Final number of systems saved into dataset: {count}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build the extended lipid form-factor/density training dataset.")
    parser.add_argument(
        "-o", "--output",
        default="lipid_dataset_extended.h5",
        help="Path to the output HDF5 file (default: lipid_dataset_extended.h5).",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Recompute and overwrite systems already present in the output file "
        "(default: skip systems already saved).",
    )
    parser.add_argument(
        "--force-id",
        action="append",
        default=[],
        help="System ID to recompute and overwrite even without --overwrite "
        "(repeatable, e.g. --force-id 771 --force-id 345).",
    )
    args = parser.parse_args()

    recompute_extended_ff_dataset(
        args.output,
        hydration_threshold=0,
        recompute_centering=False,
        small_trajs_only=False,
        overwrite=args.overwrite,
        force_ids=set(args.force_id),
    )
