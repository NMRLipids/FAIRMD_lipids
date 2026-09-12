#!/usr/bin/env python3

"""
Script for processing lipid databank systems and saving analysis results to HDF5.

This script filters systems by water-to-lipid ratio, calculates density
profiles using maicos, and exports the data for machine learning training.
"""

import logging
import os

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
from fairmd.lipids.molecules import Lipid


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
    except:
        print(f"System {system} - can't load ApL")
        no_error_flag = False
        apl = -1

    try:
        thickness = get_thickness(system)
    except:
        print(f"System {system} - can't load thickness")
        no_error_flag = False
        thickness = -1
    return apl, thickness, no_error_flag


def center_trajectory(
    system: System,
    uc: UniverseConstructor,
    last_atom,
    g3_atom,
    eq_time,
    logger,
    *,
    recompute: bool = False,
) -> str:
    """
    Centers the simulation trajectory for analysis, handling different software backends.

    Coordinates trajectory centering using either Gromacs commands or MDAnalysis
    (sequential or parallel) based on the simulation metadata and environment
    configuration.

    Args:
        u (MDAnalysis.Universe): The initial MDAnalysis universe.
        uc (UniverseConstructor): Object containing simulation path and topology info.
        spath (str): Full path to the simulation directory.
        last_atom (int/str): Index or name of the last atom for centering reference.
        g3_atom (int/str): Index or name of the glycerol 3 atom for orientation.
        eq_time (float): Equilibration time to skip in milliseconds.
        logger (logging.Logger): Logger instance for status and error reporting.

    Returns:
    str: The file path to the newly created centered trajectory file
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

    Returns:
        tuple: A triplet of space-separated strings (head_atoms, tail_atoms, backbone_atoms)
            containing the atom names for each respective fragment.
    """
    head_atoms = ""
    tail_atoms = ""
    backbone_atoms = ""
    for atom in mapping_dict:
        fragment = mapping_dict[atom]["FRAGMENT"]
        atom_name = mapping_dict[atom]["ATOMNAME"]
        if fragment == "headgroup":
            head_atoms += atom_name + " "
        elif fragment == "sn-1" or fragment == "sn-2" or fragment == "tail":
            tail_atoms += atom_name + " "
        elif fragment == "glycerol backbone":
            backbone_atoms += atom_name + " "
        else:
            print(f"Invalid atom - {atom} - {atom_name} - {fragment}")
    return (head_atoms, tail_atoms, backbone_atoms)


def create_fragment_selectors(lipid_names):
    """
    Create MDAnalysis atom selection strings for lipid fragments across a system.

    Iterates through a list of lipid names, registers their fragment mappings,
    and constructs 'name ...' strings used to select specific fragments in
    a simulation.

    Args:
        lipid_names (list of str): List of lipid molecule names present in the system.

    Returns:
        list of str: A list containing three selection strings in the order:
            [head_selector, tail_selector, backbone_selector].
    """
    head_selector, tail_selector, backbone_selector = "name ", "name ", "name "
    for lipid in lipid_names:
        lipid_class = Lipid(lipid)
        lipid_class.register_mapping()

        mapping_dict = lipid_class.mapping_dict
        head_atoms, tail_atoms, backbone_atoms = separate_lipid_atoms(mapping_dict)

        head_selector += head_atoms
        tail_selector += tail_atoms
        backbone_selector += backbone_atoms

    return [head_selector, tail_selector, backbone_selector]


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

    def save_system(self, system, scalar_data, form_factor, total_dens, mol_densities, frag_densities):
        """
        Save a single system's results to the HDF5 file.

        Organizes data into groups for metadata, axes, form factors, and
        various electron density profiles (total, fragment-based, and molecule-based).
        If the system ID already exists, it is overwritten to ensure data integrity.

        Args:
            system (dict): System metadata from the databank.
            scalar_data (dict): Dictionary containing 'ApL' and 'thickness'.
            form_factor (tuple): (q_pos, profile, dprofile) for the form factor.
            total_dens (tuple): (r_pos, profile, dprofile) for total electron density.
            mol_densities (list of tuples): List of (profile, dprofile) for each molecule type.
            frag_densities (list of tuples): List of (profile, dprofile) for lipid fragments.
        """
        with h5py.File(self.filename, "a") as f:
            sys_id = str(system["ID"])

            if sys_id in f:
                print(f"Warning: System {sys_id} already exists in HDF5. Overwriting.")
                del f[sys_id]

            grp = f.create_group(sys_id)

            grp.attrs["path"] = system.get("path", "")
            grp.attrs["ApL"] = scalar_data.get("ApL", 0)
            grp.attrs["thickness"] = scalar_data.get("thickness", 0)

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
            for i, (profile, dprofile) in enumerate(mol_densities):
                sub_grp = mol_grp.create_group(f"mol_{i}")
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
) -> None:
    """
    Recompute the extended form factor dataset for all systems in the databank.

    This function iterates through all systems, checks their suitability, and
    processes them to extract form factors and density profiles. The results
    are saved into an HDF5 file using the HDF5LipidWriter class. Systems that
    do not meet the criteria (e.g., water-to-lipid ratio) are skipped.
    """
    systems = initialize_databank()
    logger = logging.getLogger(__name__)
    writer = HDF5LipidWriter(h5fpath)

    count = 0
    print(f"Number of systems: {len(systems)}")
    for system in systems:
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
        except:
            continue

        eq_time = float(system["TIMELEFTOUT"]) * 1000
        last_atom, g3_atom = first_last_carbon(system, logger)

        u = center_trajectory(system, uc, last_atom, g3_atom, eq_time, logger, recompute=recompute_centering)

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

        molecule_types_selector = [
            f"resname {system['COMPOSITION'][molkey]['NAME']}" for molkey in system.content
        ]
        dens_molecule = []
        for selector in molecule_types_selector:
            print(f"Calculating {selector.split(' ')[1]} density")
            molecule_group = u.select_atoms(selector)
            dens_molecule_runner = DensityPlanar(
                molecule_group,
                dens="electron",
                **dens_options,
            ).run()
            dens = (dens_molecule_runner.results.profile, dens_molecule_runner.results.dprofile)
            dens_molecule.append(dens)

        fragment_selectors = create_fragment_selectors(system.lipids.keys())
        dens_fragment = []
        frag_labels = ["head", "tail", "backbone"]
        for i, selector in enumerate(fragment_selectors):
            print(f"Calculating {frag_labels[i]} density")
            fragment_group = u.select_atoms(selector)
            dens_fragment_runner = DensityPlanar(
                fragment_group,
                dens="electron",
                **dens_options,
            ).run()
            dens = (dens_fragment_runner.results.profile, dens_fragment_runner.results.dprofile)
            dens_fragment.append(dens)

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
    h5fpath = "lipid_dataset_extended.h5"
    recompute_extended_ff_dataset(h5fpath, hydration_threshold=0, recompute_centering=False, small_trajs_only=False)