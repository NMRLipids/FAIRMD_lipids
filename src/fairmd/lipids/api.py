"""
API functions for analyzing the FAIRMD Lipids Databank.

Functionality are organized into few groups:

1. Class :class:`UniverseConstructor`, which help one to create MDAnalysis.Universe from Databank's
   :class:`System <fairmd.lipids.core.System>`
2. Functions that extract computed properties:

   - :func:`get_OP`
   - :func:`get_thickness`
   - :func:`get_eqtimes`
   - :func:`get_FF`
3. Functions that extract post-processed properties:

   - :func:`get_mean_ApL`
   - :func:`get_total_area`
4. Auxiliary functions for better interface with *MDAnalysis*

   - :func:`mda_gen_selection_mols`
   - :func:`mda_read_trj_tilt_angles`
"""

import json
import logging
import math
import os
import sys
import warnings
from collections.abc import Container

import MDAnalysis as mda
import numpy as np

from fairmd.lipids import FMDL_SIMU_PATH
from fairmd.lipids.core import System
from fairmd.lipids.databankio import download_resource_from_uri, resolve_file_url
from fairmd.lipids.molecules import Molecule, lipids_set
from fairmd.lipids.schema_validation.engines import get_struc_top_traj_fnames

logger = logging.getLogger(__name__)


def get_thickness(system: System) -> float:
    """
    Get thickness for a simulation defined with ``system`` from the ``thickness.json``.

    :param system: Simulation object.

    :return: membrane thickess (nm) or raise exception
    """
    thickness_path = os.path.join(FMDL_SIMU_PATH, system["path"], "thickness.json")
    try:
        with open(thickness_path) as f:
            thickness = json.load(f)
        thickness_v = float(thickness)
    except FileNotFoundError:
        print("No thickness information for system#{}.".format(system["ID"]), file=sys.stderr)
        raise
    except ValueError:
        print("Thickness information for system#{} is invalid.".format(system["ID"]), file=sys.stderr)
        raise
    else:
        return thickness_v


def get_eqtimes(system: System) -> dict:
    """
    Return relative equilibration time for each lipid of ``system``.

    :param system: Simulation object.

    :return: dictionary of relative equilibration times for each lipid
    """
    eq_times_path = os.path.join(FMDL_SIMU_PATH, system["path"], "eq_times.json")

    try:
        with open(eq_times_path) as f:
            eq_time_dict = json.load(f)
    except FileNotFoundError as e:
        msg = "No equilibration time information for system#{}.".format(system["ID"])
        raise FileNotFoundError(msg) from e
    except json.JSONDecodeError as e:
        msg = "Equilibration times information for system#{} is invalid.".format(system["ID"])
        raise ValueError(msg) from e

    return eq_time_dict


def get_OP(system: System) -> dict:  # noqa: N802 (API name)
    """
    Return a dictionary with the order parameter data for each lipid in ``system``.

    :param system: Simulation object.

    :return: dictionary contaning, for each lipid, the order parameter data:
             average OP, standard deviation, and standard error of mean. Contains
             None if ``LipidNameOrderParameters.json`` missing.
    """
    sim_op_data = {}  # order parameter data for each type of lipid
    for mol in system["COMPOSITION"]:
        if mol not in lipids_set:
            continue
        fname = os.path.join(
            FMDL_SIMU_PATH,
            system["path"],
            mol + "OrderParameters.json",
        )
        # it always returns dictionary but values can be empty
        if not os.path.isfile(fname):
            warnings.warn(f"{fname} not found for {system['ID']}", stacklevel=2)
            sim_op_data[mol] = None
            continue
        op_data = {}
        try:
            with open(fname) as json_file:
                op_data = json.load(json_file)
        except json.JSONDecodeError as e:
            msg = f"Order parameter data in {fname} is invalid for {system['ID']}"
            raise ValueError(msg) from e
        sim_op_data[mol] = op_data
    return sim_op_data


def get_FF(system: System) -> np.ndarray:  # noqa: N802 (API name)
    """
    Get numpy table of FormFactor curve.

    :param system: Simulation object
    :return: (q,FF,err) numpy table
    """
    fn = os.path.join(
        FMDL_SIMU_PATH,
        system["path"],
        "FormFactor.json",
    )
    try:
        with open(fn) as json_file:
            sim_ff_data = json.load(json_file)
    except FileNotFoundError as e:
        msg = "The form-factor data is missing for system#{}.".format(system["ID"])
        raise FileNotFoundError(msg) from e
    except json.JSONDecodeError as e:
        msg = "The form-factor data for system#{} is invalid.".format(system["ID"])
        raise ValueError(msg) from e
    return np.array(sim_ff_data)


def get_quality(
    system: System,
    *,
    part: str = "total",
    lipid: str | None = None,
    experiment: str = "both",
) -> dict:
    """
    Return quality metrics for a given system.

    :param system: Simulation system
    :param part: Part of the system to evaluate quality for (total|tails|headgroup)".
    :param lipid: Lipid name to evaluate quality for (if None, evaluates for all lipids).
    :param experiment: Experiment type to evaluate quality against (XR|NMR|both).
    """
    if part not in ["total", "headgroup", "tails"] or experiment not in ["FF", "OP", "both"]:
        msg = "Invalid values for `part` or `experiment`!"
        raise ValueError(msg)
    if (part != "total" or lipid is not None) and experiment in ["both", "FF"]:
        msg = "Combined or form-factor qualities are available only for the entire system!"
        raise ValueError(msg)
    if part == "total" and experiment == "both":
        msg = "Quality for both experiments is not implemented!"
        raise NotImplementedError(msg)

    q = -100

    if experiment == "FF":
        spath = os.path.join(FMDL_SIMU_PATH, system["path"], "FormFactorQuality.json")
        if not os.path.isfile(spath):
            return np.nan
        with open(spath) as fd:
            q = float(json.load(fd)[0])

    if experiment == "OP" and lipid is None:
        spath = os.path.join(FMDL_SIMU_PATH, system["path"], "SYSTEM_quality.json")
        if not os.path.isfile(spath):
            return np.nan
        with open(spath) as fd:
            qdict = json.load(fd)
        return float(qdict[part])

    if experiment == "OP" and lipid is not None:
        if lipid not in system.content:
            msg = f"Lipid {lipid} is not in the system composition!"
            raise ValueError(msg)
        spath = os.path.join(FMDL_SIMU_PATH, system["path"], lipid + "_FragmentQuality.json")
        if not os.path.isfile(spath):
            return np.nan
        with open(spath) as fd:
            qdict = json.load(fd)
        if part == "tails":
            vals = [qdict.get("sn-1"), qdict.get("sn-2")]
            vals = [v for v in vals if v is not None]

            if not vals:
                return np.nan

            return float(np.nanmean(vals))

        else:
            return float(qdict[part])
    return q


def get_mean_ApL(system: System) -> float:  # noqa: N802 (API name)
    """
    Calculate average area per lipid for a system.

    :param system: Simulation object.

    :return: area per lipid (Å^2)
    """
    path = os.path.join(FMDL_SIMU_PATH, system["path"], "apl.json")
    try:
        with open(path) as f:
            data = json.load(f)
    except FileNotFoundError as e:
        msg = "Area per lipid data is absent for system #{}".format(system["ID"])
        raise FileNotFoundError(msg) from e
    except json.JSONDecodeError as e:
        msg = "Area per lipid data for system #{} in {} is invalid.".format(system["ID"], path)
        raise ValueError(msg) from e
    vals = np.array(list(data.values()))
    return vals.mean()


def get_total_area(system: System) -> float:
    """
    Return area of the membrane in the simulation box.

    :param system: Simulation object.

    :return: area of the system (Å^2)
    """
    apl = get_mean_ApL(system)
    return system.n_lipids * apl / 2


class UniverseConstructError(Exception):
    """Specific error for UniverseConstructor"""


class UniverseConstructor:
    """
    Class operating with downloading and constructing Universe for the :class:`System <fairmd.lipids.core.System>`.

    To use this class, one instantinate it with a particular system, and then download.

    .. code-block:: python

        s = systems.loc(120)
        uc = UniverseConstructor(s)
        uc.download_mddata()

    After this, the pointer :attr:`uc.paths <paths>` will show which files are available to work with.
    Finally, you can run :meth:`uc.build_universe() <build_universe>` to get the ``MDAnalysis.Universe`` object.
    """

    def __init__(self, s: System) -> None:
        """
        Create an empty instance. No auto-download or auto-universe.

        :param s: Simulation object.
        """
        self._s = s
        self._paths = {
            "struc": None,
            "top": None,
            "traj": None,
            "energy": None,
        }

    @property
    def system(self) -> System:
        """Link to simulation object."""
        return self._s

    @property
    def paths(self) -> dict[str, str | None]:
        """Return dicts of absolute paths of downloaded files.

        Allowed fields are: *struc*, *traj*, *top*, *energy*. If they are not ``None``, then
        the corresponding file is downloaded and the full path is the value.
        """
        return self._paths

    def download_mddata(self, *, skip_traj: bool = False) -> None:
        """
        Download all the files. Previously downloaded are skipped.

        :param skip_traj: Download only TOP&struc for further constructing single-frame universe
        """
        gpath = os.path.join(FMDL_SIMU_PATH, self._s["path"])
        struc, top, trj = get_struc_top_traj_fnames(self._s)

        def _resolve_dwnld(fname: str) -> str:
            fpath = os.path.join(gpath, fname)
            if self._s["DOI"] == "localhost":
                if not os.path.isfile(fpath):
                    msg = f"File {fpath} must be predownloaded for {self._s}"
                    raise FileNotFoundError(msg)
                return fpath
            if os.path.isfile(fpath):
                # do not download if exists
                return fpath
            url = resolve_file_url(self._s["DOI"], fname)
            _ = download_resource_from_uri(url, fpath, max_restarts=5)
            return fpath

        if struc is not None:
            self._paths["struc"] = _resolve_dwnld(struc)
        if top is not None:
            self._paths["top"] = _resolve_dwnld(top)
        if trj is not None and not skip_traj:
            self._paths["traj"] = _resolve_dwnld(trj)

    def clear_mddata(self) -> None:
        """Clear downloaded MD data. For DOI=localhost, do nothing."""
        if self._s["DOI"] == "localhost":
            for k in self._paths:
                self._paths[k] = None
            return
        for k, v in self._paths.items():
            if v is None:
                continue
            print(f"Clearing {k}-file..", end="", flush=True)
            os.remove(v)
            print("OK")
            self._paths[k] = None

    def build_universe(self) -> mda.Universe:
        """Build MDAnalysis Universe.

        Replaces outdated `system2MDanalysisUniverse`.
        """
        if not any(self._paths.values()):
            msg = "You **MUST** run `download_mddata` before `build_universe`"
            raise UniverseConstructError(msg)

        if self._paths["top"] is None:
            u = (
                mda.Universe(self._paths["struc"])
                if self._paths["traj"] is None
                else mda.Universe(self._paths["struc"], self._paths["traj"])
            )
        else:
            try:
                if self._paths["traj"] is None:
                    u = (
                        mda.Universe(self._paths["top"])
                        if self._paths["struc"] is None
                        else mda.Universe(self._paths["top"], self._paths["struc"])
                    )
                else:
                    u = mda.Universe(self._paths["top"], self._paths["traj"])
            except OSError as e:  # exception for corrupted topology file
                print(
                    f"We got exception.. == \n{e}\n == ..and assume that TOPOLOGY is file is corrupted", file=sys.stderr
                )
                if self._paths["struc"] is None:
                    msg = "TOPOLOGY is corrupted, and no STRUCTURE is given"
                    raise UniverseConstructError(msg) from e
                u = (
                    mda.Universe(self._paths["struc"])
                    if self._paths["traj"] is None
                    else mda.Universe(self._paths["struc"], self._paths["traj"])
                )
            # other exceptions are raised as is
        return u


def mda_gen_selection_mols(system: System, molecules: Container[Molecule] | None = None) -> str:
    """
    Return a MDAnalysis selection string covering all the molecules (default None means "lipids").

    :param system: FAIRMD Lipids dictionary defining a simulation.
    :param molecules: container of molecule objects to be included in the selection.

    :return: a string using MDAnalysis notation that can used to select all lipids from
             the ``system``.
    """
    res_set = set()
    molecules = system.lipids.values() if molecules is None else molecules
    for key, mol in system.content.items():
        if mol in molecules:
            try:
                for atom in mol.mapping_dict:
                    res_set.add(mol.mapping_dict[atom]["RESIDUE"])
            except (KeyError, TypeError):
                res_set.add(system["COMPOSITION"][key]["NAME"])
    sorted_res = sorted(res_set)
    return "resname " + " or resname ".join(sorted_res)


def mda_read_trj_tilt_angles(
    resname: str,
    a_name: str,
    b_name: str,
    universe: mda.Universe,
):
    """
    Calculate the AB vector angles with respect to membrane normal from the simulation.

    :param resname: residue name of the molecule for which the P-N vector angle will be calculated
    :param a_name: name of the A atom in the simulation
    :param b_name: name of the B atom in the simulation
    :param universe: MDAnalysis universe of the simulation to be analyzed

    :return: tuple (angles of all molecules as a function of time,
                    time averages for each molecule,
                    the average angle over time and molecules,
                    the error of the mean calculated over molecules)
    """

    # Auxiliary internal function
    def _calc_angle(atoms, com) -> float:
        """
        :meta private:
        Calculate the angle between the vector and z-axis in degrees.

        No PBC check! Calculates the center of mass of the selected atoms to invert bottom leaflet vector
        """
        vec = atoms[1].position - atoms[0].position
        d = math.sqrt(np.square(vec).sum())
        cos = vec[2] / d
        # values for the bottom leaflet are inverted so that
        # they have the same nomenclature as the top leaflet
        cos *= math.copysign(1.0, atoms[0].position[2] - com)
        try:
            angle = math.degrees(math.acos(cos))
        except ValueError:
            if abs(cos) >= 1.0:
                print(f"Cosine is too large = {cos} --> truncating it to +/-1.0")
                cos = math.copysign(1.0, cos)
                angle = math.degrees(math.acos(cos))
        return angle

    selection = universe.select_atoms(
        "resname " + resname + " and (name " + a_name + ")",
        "resname " + resname + " and (name " + b_name + ")",
    ).atoms.split("residue")
    com = universe.select_atoms(
        "resname " + resname + " and (name " + a_name + " or name " + b_name + ")",
    ).center_of_mass()

    n_res = len(selection)
    n_frames = len(universe.trajectory)
    angles = np.zeros((n_res, n_frames))

    res_aver_angles = [0] * n_res
    res_std_error = [0] * n_res
    j = 0

    for _ in universe.trajectory:
        for i in range(n_res):
            residue = selection[i]
            angles[i, j] = _calc_angle(residue, com[2])
        j = j + 1
    for i in range(n_res):
        res_aver_angles[i] = sum(angles[i, :]) / n_frames
        res_std_error[i] = np.std(angles[i, :])

    total_average = sum(res_aver_angles) / n_res
    total_std_error = np.std(res_aver_angles) / np.sqrt(n_res)

    return angles, res_aver_angles, total_average, total_std_error


# -------------------------------------- SEPARATED PART (??) ----------------------
