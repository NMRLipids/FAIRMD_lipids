"""
Module for calculation of order parameters of lipid bilayers from a MD trajectory

**Authors:**
- Made by Joe,  Last edit 2017/02/02
- Refactored for performance by Gemini
"""

import re
import warnings

import MDAnalysis as mda
import numpy as np

from fairmd.lipids import progress

# Maximum bond length for a C-H bond to be considered reasonable.
bond_len_max = 1.5  # in Angstrom
bond_len_max_sq = bond_len_max**2


class _OrderParameter:
    """
    Atomic dipole with order parameter.

    Allows to store and manipulating OP related metadata (definition, name, etc.), OP trajectories,
    and methods to evaluate OPs.
    """

    def __init__(
        self,
        resname: str,
        atom_name_a: str,
        atom_name_b: str,
        univ_atom_name_a: str,
        univ_atom_name_b: str,
    ) -> None:
        """Initialize the OrderParameter object.

        It doesn't matter which atom (A or B) comes first for the OP calculation.

        :param resname: Name of the residue the atoms are in.
        :param atom_name_a: Name of the first atom in the topology.
        :param atom_name_b: Name of the second atom in the topology.
        :param univ_atom_name_a: Generic/mapping name for atom A.
        :param univ_atom_name_b: Generic/mapping name for atom B.
        :raises RuntimeError: If any of the provided names are empty strings.
        """
        self.resname = resname
        self.aname_a = atom_name_a
        self.aname_b = atom_name_b
        self.m_aname_a = univ_atom_name_a
        self.m_aname_b = univ_atom_name_b
        self.name = f"{univ_atom_name_a} {univ_atom_name_b}"

        for field_name, field_value in self.__dict__.items():
            if isinstance(field_value, str):
                if not field_value.strip():
                    msg = (
                        f"Provided name for field '{field_name}' is empty! "
                        "Cannot use empty names for atoms and OP definitions."
                    )
                    raise RuntimeError(msg)
            else:
                msg = f"Provided value for '{field_name}' is not a string: {field_value}. "
                raise TypeError(msg)

        self._avg = None
        self._std = None
        self._stem = None

        self.traj = []  # For storing final OP results.
        self.selection = []  # List of AtomGroups, one for each residue.
        self.atomgroup = None  # A single AtomGroup containing all atoms for this OP.

    def finalize(self) -> None:
        """Finalize the OP object by calculating average, stddev, and stem."""
        n = len(self.traj)
        if n == 0:
            msg = f"No trajectory data available for OP '{self.name}'. Cannot finalize."
            raise RuntimeError(msg)
        self._std = np.std(self.traj)
        self._avg = np.mean(self.traj)
        self._stem = self._std / np.sqrt(n - 1) if n > 1 else 0


    @property
    def avg_std_stem(self) -> tuple[float, float, float]:
        """Average, stddev, and standard error of the mean of OPs."""
        return self._avg, self._std, self._stem


def _read_trajs_calc_OPs(  # noqa: N802
    op_obj_list: list[_OrderParameter],
    universe: mda.Universe,
) -> None:
    """Read a Universe trajectory and calculate Order Parameters ("S").

    This function calculates the order parameters for each definition in ``op_obj_list``.
    This version is optimized for single-core performance using vectorized calculations.
    The results are stored in-place in the ``traj`` attribute of the objects in ``op_obj_list``.

    :param op_obj_list: A list of _OrderParameter objects to be processed.
    :param universe: MDAnalysis Universe containing topology and trajectory.
    """
    # --- 1. Setup Universe and Atom Selections ---
    mol = universe
    improper_ops = []

    for i, op in enumerate(op_obj_list):
        sel_str = f"resname {op.resname} and name {op.aname_a} {op.aname_b}"
        selection_by_residue = mol.select_atoms(sel_str).split("residue")

        if not selection_by_residue:
            warnings.warn(
                f"Selection is empty: [{sel_str}]. Check residue and atom names in the mapping file.",
                UserWarning,
                stacklevel=2,
            )
            improper_ops.append(i)
            continue

        # Validate that each residue selection contains exactly two atoms
        valid_selection = []
        for res in selection_by_residue:
            if res.n_atoms != 2:  # noqa: PLR2004
                warnings.warn(
                    f"Selection 'name {op.aname_a} {op.aname_b}' in residue "
                    f"{res.resids[0]} contains {res.n_atoms} atoms, but should be 2. "
                    "This residue will be skipped.",
                    UserWarning,
                    stacklevel=2,
                )
            else:
                valid_selection.append(res)

        if not valid_selection:
            warnings.warn(
                f"No valid atom pairs found for selection: [{sel_str}]",
                UserWarning,
                stacklevel=2,
            )
            improper_ops.append(i)
            continue

        op.selection = valid_selection
        # Create a single, combined AtomGroup for efficient, vectorized access
        op.atomgroup = mda.AtomGroup([atom for res in op.selection for atom in res])

    # Remove OP definitions that resulted in invalid selections
    for i in sorted(improper_ops, reverse=True):
        del op_obj_list[i]

    n_frames = len(mol.trajectory)
    for op in op_obj_list:
        n_res = len(op.selection)
        # We accumulate sums here, so initialize with zeros
        op.traj = np.zeros(n_res, dtype=np.float64)

    print("Processing trajectory with optimized single-core engine...")
    for _ in progress(
        mol.trajectory,
        total=n_frames,
        unit="frame",
        desc="Processing trajectory",
    ):
        for op in op_obj_list:
            if op.atomgroup is None or len(op.atomgroup) == 0:
                continue

            # Get all atom positions for this OP in one go
            # Shape is n_residues*2 x 3
            positions = op.atomgroup.positions

            # Reshape to easily access atom pairs
            # Shape: (n_residues, 2, 3) where axis 1 is [atom_A, atom_B]
            positions_reshaped = positions.reshape(-1, 2, 3)

            # Calculate vectors between atom pairs for all residues at once
            vec = positions_reshaped[:, 1, :] - positions_reshaped[:, 0, :]

            # Calculate squared distance for all residues
            d2 = np.sum(vec**2, axis=1)

            # Create a mask for valid bond lengths to avoid unnecessary calculations
            # and warnings for atoms that are too far apart (e.g., due to PBC issues).
            valid_mask = d2 <= bond_len_max_sq

            # Initialize cos2 array. Invalid long bonds remain zero and are
            # excluded by the existing valid-mask policy.
            cos2 = np.zeros_like(d2)

            # Safely calculate cosine-squared of the angle with the z-axis
            # for all valid vectors simultaneously.
            # Zero-length pairs are represented as NaN: their direction is
            # undefined and must not silently contribute to the result.
            d2_valid = d2[valid_mask]
            vec_valid = vec[valid_mask]
            cos2[valid_mask] = np.divide(
                vec_valid[:, 2] ** 2,
                d2_valid,
                out=np.full_like(d2_valid, np.nan),
                where=d2_valid != 0,
            )

            # Calculate order parameters for all residues
            op_values = 0.5 * (3.0 * cos2 - 1.0)

            # Add the results for the current frame to the running sum.
            # Invalid long bonds remain 0; undefined zero-length pairs make
            # the accumulated result NaN, exposing the bad input.
            op.traj += op_values

    # Average the accumulated sums over all frames
    for op in op_obj_list:
        if n_frames > 0:
            op.traj /= n_frames
        # Convert back to a list to maintain original API behavior
        op.traj = op.traj.tolist()
        op.finalize()


def _parse_op_input(mapping_dict: dict, lipid_resname: str) -> list[_OrderParameter]:
    """Parse a mapping dictionary to form a list of C-H pairs for OP calculation.

    :param mapping_dict: The mapping dictionary.
    :param lipid_resname: The default lipid residue name.
    :return: A list of _OrderParameter instances.
    """
    opvals = []
    atom_c = []
    atom_h = []
    resname = lipid_resname

    # Regex to identify carbon and hydrogen atoms from mapping keys
    regexp_h = re.compile(r"M_([A-Z0-9]*C[0-9]*|G[0-9]*|C[0-9]*)H[0-9]*_M")
    regexp_c = re.compile(r"M_([A-Z0-9]*C[0-9]*|G[0-9]{1,2}|C[0-9]{1,2})_M")

    for mapping_key, value in mapping_dict.items():
        if not isinstance(value, dict) or "ATOMNAME" not in value:
            continue

        if regexp_c.search(mapping_key) and not regexp_h.search(mapping_key):
            atom_c = [mapping_key, value["ATOMNAME"]]
            # Use residue name from mapping if available, otherwise use default
            resname = value.get("RESIDUE", lipid_resname)
            atom_h = []  # Reset hydrogen atom
        elif regexp_h.search(mapping_key):
            atom_h = [mapping_key, value["ATOMNAME"]]
        else:
            atom_c, atom_h = [], []  # Reset for non-matching keys

        if atom_h and not atom_c:
            warnings.warn(
                f"Cannot define carbon for the hydrogen {atom_h[0]} ({atom_h[1]})",
                UserWarning,
                stacklevel=2,
            )
            continue

        # If both a carbon and a hydrogen have been found, create the pair
        if atom_h and atom_c:
            op = _OrderParameter(resname, atom_c[1], atom_h[1], atom_c[0], atom_h[0])
            opvals.append(op)
            # Important: Reset hydrogen to look for the next one for the same carbon
            atom_h = []

    return opvals


def find_OP(  # noqa: N802
    mdict: dict,
    universe: mda.Universe,
    lipid_name: str,
) -> list[_OrderParameter]:
    """Externally used function for computing OP values.

    :param mdict: The mapping dictionary.
    :param universe: MDAnalysis Universe containing topology and trajectory.
    :param lipid_name: The residue name of the lipid.

    :return: A list of _OrderParameter instances with calculated data.
    """
    op_pairs = _parse_op_input(mdict, lipid_name)
    _read_trajs_calc_OPs(op_pairs, universe)
    return op_pairs
