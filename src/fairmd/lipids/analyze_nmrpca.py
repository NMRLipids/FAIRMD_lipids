"""
NMRPCA module calculate the relaxation time based on the PCAlipids analysis.

For details check:

Principal Component Analysis of Lipid Molecule Conformational Changes in Molecular Dynamics Simulations
Pavel Buslaev, Valentin Gordeliy, Sergei Grudinin, and Ivan Gushchin.
*Journal of Chemical Theory and Computation* (2016), 12(3), 1019–1028.
DOI: `10.1021/acs.jctc.5b01106 <https://doi.org/10.1021/acs.jctc.5b01106>`_

and

Principal component analysis highlights the influence of temperature,
curvature and cholesterol on conformational dynamics of lipids
Pavel Buslaev, Khalid Mustafin, and Ivan Gushchin
*Biochimica et Biophysica Acta (BBA) - Biomembranes*
Volume 1862, Issue 7, 1 July 2020, 183253
DOI: `10.1016/j.bbamem.2020.183253 <https://doi.org/10.1016/j.bbamem.2020.183253>`_

The main idea is to run PCA on concatenated lipid trajectory,
calculate the autocorrelation times, which are linearly related
to equilibration times of individual lipids. The parameters of linear
transform are calculated based on the trajectories from NMRlipids
databank.

The initial code was developed by
Dr. Pavel Buslaev (pavel.i.buslaev@jyu.fi),
Dr. Samuli Olilla,
Patrik Kula,
Alexander Kuzmin

.. envvar:: MEM_CUTOFF
   :noindex:

   Environmental variable which sets the trajectory size cutoff (in Bytes)
   after which the algorithm will start skipping frames for in-memory
   representation for computing ACF.

   Default: ``1,500,000,000``.
"""

import json
import os
import warnings

import MDAnalysis as mda
import numpy as np
from deprecated import deprecated
from MDAnalysis.analysis.base import AnalysisFromFunction
from scipy import signal

from fairmd.lipids import FMDL_SIMU_PATH
from fairmd.lipids.core import System
from fairmd.lipids.databankio import download_resource_from_uri, resolve_download_file_url
from fairmd.lipids.databankLibrary import lipids_set
from fairmd.lipids.settings.molecules import Lipid

# suppress some MDAnalysis warnings issued from mda.analysis.align
warnings.filterwarnings("ignore")
from MDAnalysis.analysis import align  # noqa: E402

# TODO: now there are only regular phospholipids.
# The list should be verified by method authors.
ALLOWLIPIDS = [
    "POPC",
    "POPG",
    "POPS",
    "POPE",
    "PYPC",
    "DMPC",
    "DPPC",
    "DPPE",
    "DPPG",
    "DEPC",
    "DRPC",
    "DYPC",
    "DLPC",
    "DLIPC",
    "DOPC",
    "DOPE",
    "DDOPC",
    "DOPS",
    "DSPC",
    "DAPC",
    "SDPE",
    "SOPC",
    "POPI",
    "SAPI",
    "SAPI24",
    "SAPI25",
    "SLPI",
]

merge_cutoff = 2.0
trj_size_cutoff = (
    1500000000  # This is the boundary set empirically by @pbuslaev
    if "MEM_CUTOFF" not in os.environ
    else int(os.environ["MEM_CUTOFF"])
)

TAILSN1 = "sn-1"
TAILSN2 = "sn-2"
HEADGRP = "headgroup"

# In case you want to test for a specific coordinate change the flag to yes
TEST = False


class Parser:
    """
    Class Parser is a basic class to work with the trajectory. It has basic utility
    for preparing trajectories:

    1. Checking if simulation has a correct name
    2. Checking if trajectory is downloaded. Downloading, if not.
    3. Calling AmberMerger to merge head groups and tails for Amber
       trajectories
    4. Concatenate trajectories

    :param system: simulation data
    :param eq_time_fname: name of the output file
    :param path: name of a particular trajectory
    :param v: verbosity for testing
    """

    def __init__(self, system: System, eq_time_fname="eq_times.json", path=None, v=True):
        # Technical
        self.verbose = v
        self.error = 0

        # Path
        self.root = FMDL_SIMU_PATH
        self.eq_time_fname = eq_time_fname

        # Extracting data from readme
        self._path = os.path.join(self.root, system["path"])
        print("Indexing path:", self._path)
        if self.verbose:
            print(f"Parser: Processing trajectory {self._path}")
        self.doi = system["DOI"]
        self.soft = system["SOFTWARE"]
        if self.soft == "openMM" or self.soft == "NAMD":
            try:
                self.trj = system["TRJ"][0][0]
                self.tpr = system["PDB"][0][0]
            except Exception:
                print("Parser: Did not find trajectory or pdb for openMM or NAMD")
                self.error = 1
        else:
            try:
                self.trj = system["TRJ"][0][0]
                self.tpr = system["TPR"][0][0]
            except Exception:
                print("Parser: Did not find trajectory or tpr for GROMACS")
                self.trj = ""
                self.tpr = ""
                self.error = 1
        self.trj_name = os.path.join(self._path, self.trj)
        self.tpr_name = os.path.join(self._path, self.tpr)
        self.trj_len = system["TRJLENGTH"] / 1000  # ns

        self.size = system["TRAJECTORY_SIZE"]

        self.composition = system["COMPOSITION"]
        self.lipids: dict[str, Lipid] = {k: v for k, v in system.content.items() if k in lipids_set}
        self.path = path

    @deprecated(reason="Downloading should happen somewhere else.Not in the analysis classes.")
    def validate_path(self):
        """
        Path validation. Behaviour depends on the input.

        1. If ther were any errors, validation fails. Currently the only
           error tested is that .xtc and .tpr files are not present. This
           is the case for non-GROMACS trajectories.
        2. If path is not provided, parser is iterating over all trajectories.
        3. If path is provided and current path is equal to the provided one,
           parser reports that it finds the trajectory for the analysis.

        TODO: must be removed. Substitute path-validation part.
        """
        if self.error > 0:
            # This is an error message. Printing even in silent mode
            print(f"Parser: Can't read TPR/PDB and TRJ from README for {self._path}")
        if self.verbose and not self.path:
            print(
                f"Parser: Iterating over all trajectories. Current trajectory is {self._path}",
            )
            # return
        # if self.path == self.indexingPath:
        if os.path.isfile(os.path.join(self._path, self.eq_time_fname)):
            if self.verbose:
                print("Parser: Found file with equilibration data. \nNot processing the trajectory")
                return -1
        if self.verbose:
            print(f"Parser: Found trajectory {self._path}")
        return 0

    @deprecated(reason="Downloading should happen somewhere else.Not in the analysis classes.")
    def download_traj(self):
        """
        TODO: must be removed. We don't download trajectory in the analysis part

        Basic trajectory and TPR download.
        """
        print("Downloading")
        if not os.path.isfile(self.tpr_name):
            self.tpr_url = resolve_download_file_url(self.doi, self.tpr)
            # This is a log message. Printing even in silent mode
            print("Parser: Downloading tpr ", self.doi)
            # urllib.request.urlretrieve(self.tpr_url, self.tpr_name)
            download_resource_from_uri(self.tpr_url, self.tpr_name)
        if not os.path.isfile(self.trj_name):
            self.trj_url = resolve_download_file_url(self.doi, self.trj)
            # This is a log message. Printing even in silent mode
            print("Parser: Downloading trj ", self.doi)
            # urllib.request.urlretrieve(self.trj_url, self.trj_name)
            download_resource_from_uri(self.trj_url, self.trj_name)

    def prepare_gmx_traj(self):
        """
        Preparing trajectory. If centered trajectory is found, use it. If whole
        trajectory is found, use it. Otherwise, call gmx trjconv to make whole
        trajectory. The selected trajectory is loaded into Universe.
        """
        # Look for centered.xtc
        if os.path.isfile(os.path.join(self._path, "centered.xtc")):
            print("Parser: Founder centered trajectory: centered.xtc")
            self.trj_name = os.path.join(self._path, "centered.xtc")
        # Look for whole.xtc
        elif os.path.isfile(os.path.join(self._path, "whole.xtc")):
            print("Parser: Founder whole trajectory: whole.xtc")
            self.trj_name = os.path.join(self._path, "whole.xtc")
        # Run gmx trjconv
        else:
            print("Parser: Making trajectory whole")
            trj_out_name = os.path.join(self._path, "whole.xtc")
            os.system(
                f"echo System | gmx trjconv -f {self.trj_name} -s {self.tpr_name} -pbc mol -o {trj_out_name}",
            )
            self.trj_name = trj_out_name

        # Skip frames for large trajectories
        if self.size > trj_size_cutoff:
            import math

            _skip = math.ceil(1.1 * self.size / trj_size_cutoff)
            trj_out_name = os.path.join(self._path, "short.xtc")
            os.system(
                f"echo System | gmx trjconv -f {self.trj_name} -s {self.tpr_name} -pbc mol -o {trj_out_name} -skip 10",
            )
            self.trj_name = trj_out_name

        # Create .gro file
        self.gro_name = os.path.join(self._path, "conf.gro")
        if os.path.isfile(self.gro_name):
            print("Parser: Found gro file")
        else:
            print("Parser: Making a gro file from the first frame")
            os.system(
                f"echo System | gmx trjconv -s {self.tpr_name}  -f {self.trj_name} -dump 0 -o {self.gro_name}",
            )

        # Loading trajectory
        try:
            self.traj = mda.Universe(self.tpr_name, self.trj_name)
        except Exception:
            self.traj = mda.Universe(self.gro_name, self.trj_name)

    # TODO: not tested!!
    def prepare_OpenMM_traj(self):  # noqa: N802
        print("openMM or NAMD")
        if self.size > trj_size_cutoff:
            import math

            _skip = math.ceil(1.1 * self.size / trj_size_cutoff)
            trj_out_name = os.path.join(self._path, "short.xtc")
            if os.path.isfile(trj_out_name):
                print("Parser: Short trajectory is found")
            else:
                u = mda.Universe(self.tpr_name, self.trj_name)

                with mda.Writer(trj_out_name, u.select_atoms("all").n_atoms) as mdaw:
                    for _ in u.trajectory[::10]:
                        mdaw.write(u.select_atoms("all"))

            self.trj_name = trj_out_name

        self.traj = mda.Universe(self.tpr_name, self.trj_name)

    def prepare_traj(self):
        if self.soft == "openMM" or self.soft == "NAMD":
            self.prepare_OpenMM_traj()
        else:
            self.prepare_gmx_traj()

    def concatenate_traj(self):
        """
        Create Concatenator and corresponding concatenated trajectories for
        all lipids available for the trajectory.
        """
        self.concatenated_trajs = {}
        for lipid_name, lipid in self.lipids.items():
            if lipid_name not in ALLOWLIPIDS:
                # We do not treat cholesterols
                continue
            topology = Topology(
                self.traj,
                self.composition[lipid_name]["NAME"],
                lipid.mapping_dict,
            )
            concatenator = Concatenator(topology, self.traj, self.composition[lipid_name]["NAME"])
            self.concatenated_trajs[lipid_name] = concatenator.concatenate()

    def dump_data(self, data):
        """Write data to json file."""
        with open(os.path.join(self._path, self.eq_time_fname), "w") as f:
            json.dump(data, f)


class Topology:
    """
    Class Topology is a class needed to extract lipid specific data and also to
    merge lipids from Amber trajectories, where lipid is often represented as 3
    residue types. It has the following methods:

    1. Simple constructor, which sets the force field,residue names, trajectory,
       and loads mapping_file
    2. Mapping file loader
    3. Function that outputs atomNames
    4. Checker if the merge is needed. Currently defunct
    5. Runner, that returns lists for head, tail1, tail2 orders for merging

    Currently defunct.

    :param traj: MDAnalysis trajectory
    :param lipid_resname: resname of the lipid (str)
    :param mapping_dict: preloaded mapping_dict
    """

    def __init__(self, traj, lipid_resname: str, mapping_dict: dict):
        self.lipid_resname = lipid_resname
        self.traj = traj
        self.mapping = mapping_dict

    def atom_names(self):
        """
        Extract all names of heavy atoms from the mapping
        """
        atoms = []
        for key in self.mapping:
            atom = self.mapping[key]["ATOMNAME"]
            if atom[0] == "H":
                continue
            atoms.append(atom)
        return atoms

    def is_merge_needed(self):
        """
        Checker for merge. Currently it checks if the RESIDUE key is in the
        mapping file.

        NOTE: This has to be changed if the content of mapping file changes
        """
        if "RESIDUE" not in self.mapping[list(self.mapping.keys())[0]].keys():
            return False
        resnames = []
        for key in self.mapping:
            atom = self.mapping[key]["ATOMNAME"]
            anames = self.atom_names()
            if atom in anames:
                resnames.append(self.mapping[key]["RESIDUE"])
        resnames = set(resnames)
        if len(resnames) > 1:
            return resnames
        return False

    def get_lipid_resnames(self):
        """
        Helper function that gets the residue names for lipid, if merge is needed
        """
        resnames = self.is_merge_needed()
        if self.is_merge_needed():
            return resnames
        # How can we end up here?
        # Since currently we only call this function when Merge is needed,
        # this is a place for checking if everything is ok and we can raise
        return self.lipid_resname

    def assign_resnames(self, resnames):
        """
        Helper function that finds head, sn-1 and sn-2 tails

        NOTE: currently there are only lipids with 2 different tails available in
        databank: e.g. POPC or POPG. This leads to different names of tails. This
        won't be the case for DOPC. Currently the algorithm is using that all three
        groups differ
        """
        if self.is_merge_needed():
            resname_dict = {}
            # First find headgroup
            for key in self.mapping:
                resname = self.mapping[key]["RESIDUE"]
                if self.mapping[key]["FRAGMENT"] == HEADGRP:
                    resname_dict[HEADGRP] = resname
                    break
            for key in self.mapping:
                resname = self.mapping[key]["RESIDUE"]
                if (
                    TAILSN1 not in resname_dict
                    and self.mapping[key]["FRAGMENT"] == TAILSN1
                    and not resname_dict[HEADGRP] == resname
                ):
                    resname_dict[TAILSN1] = resname
                if (
                    TAILSN2 not in resname_dict
                    and self.mapping[key]["FRAGMENT"] == TAILSN2
                    and not resname_dict[HEADGRP] == resname
                ):
                    resname_dict[TAILSN2] = resname
                if TAILSN1 in resname_dict and TAILSN2 in resname_dict:
                    break
            # TODO: add check that all resnames from input are in the dict
            return resname_dict
        # How can we end up here?
        return None

    def run_merger(self):
        """
        Find lipid tails that correspond to a particular head-group.

        NOTE: currently there are only lipids with 2 different tails available in
        databank: e.g. POPC or POPG. This leads to different names of tails. This
        won't be the case for DOPC. Currently the algorithm is using that all three
        groups differ

        TODO: get the correspondence from structure
        """
        resnames = self.get_lipid_resnames()
        resname_dict = self.assign_resnames(resnames)
        head_residues = [
            r.atoms.select_atoms("not name H* and prop mass > 0.8")
            for r in self.traj.select_atoms(
                f"not name H* and resname {resname_dict[HEADGRP]}",
            ).residues
        ]
        sn_1_residues = [
            r.atoms.select_atoms("not name H* and prop mass > 0.8")
            for r in self.traj.select_atoms(
                f"not name H* and resname {resname_dict[TAILSN1]} and "
                f"around {merge_cutoff} (resname {resname_dict[HEADGRP]} "
                "and not name H*)",
            ).residues
        ]
        sn_2_residues = [
            r.atoms.select_atoms("not name H* and prop mass > 0.8")
            for r in self.traj.select_atoms(
                f"not name H* and resname {resname_dict[TAILSN2]} and "
                f"around {merge_cutoff} (resname {resname_dict[HEADGRP]} "
                "and not name H*)",
            ).residues
        ]
        return head_residues, sn_1_residues, sn_2_residues


class Concatenator:
    """
    Class Concatenator is a class needed to concatenate trajectory for lipid types.
    It has the following methods:

    1. Simple constructor, which sets the topology,residue names, and
       trajectory
    2. :meth:`concatenateTraj` method to do the basic by lipid concatenation
    3. :meth:`alignTraj` method to perform the alignment of the concatenated trajectory
    4. The enveloping concatenate method

    :param topology: topology for lipid
    :param traj: MDAnalysis trajectory
    :param lipid_resname: lipid resname in the trajectory
    """

    def __init__(self, topology: Topology, traj, lipid_resname: str):
        self.topology = topology
        self.traj = traj
        self.lipid_resname = lipid_resname

        if self.topology.is_merge_needed():
            self.headlist, self.tail1list, self.tail2list = self.topology.run_merger()
            # TODO: raise if length of 3 lists is not equal
        else:
            self.headlist = None
            self.tail1list = None
            self.tail2list = None

    def concatenate_traj(self):
        """
        ConcatenateTraj performs basic trajectory concatination. First, it extracts
        coordinates from trajectory, next, it reshapes the coordinate array, swaps
        time and resid axes to obtain continuous trajectories of individual lipids
        (this is needed for autocorrelation time analysis), and finally merges
        individual lipid trajectories.
        """
        traj = self.traj.trajectory
        n_frames = len(traj)

        heavy_atoms_topology = self.traj.select_atoms(
            f"resname {self.lipid_resname} and not name H*",
        )

        n_atoms_lipid = len(self.topology.atom_names())

        n_lipid = heavy_atoms_topology.n_residues
        if n_lipid * n_atoms_lipid != heavy_atoms_topology.n_atoms:
            n_lipid = heavy_atoms_topology.n_atoms // n_atoms_lipid

        # Get all coordinates n_frames, n_lipid * n_atoms_lipid
        coords = AnalysisFromFunction(lambda ag: ag.positions.copy(), heavy_atoms_topology).run().results["timeseries"]
        # Adding axis to separate individual lipid trajectories
        coords = coords.reshape((n_frames, n_lipid, n_atoms_lipid, 3))
        # Swapping time frame with lipid axis
        # Now we have continuous lipid trajectory
        coords = np.swapaxes(coords, 0, 1)
        # Reshaping to desired shape
        coords = coords.reshape((n_frames * n_lipid, n_atoms_lipid, 3))

        # Creating new Universe from the concatenated coordinates
        atom_resindex = np.repeat(0, n_atoms_lipid)
        concatenated_traj = mda.Universe.empty(
            n_atoms_lipid,
            1,
            atom_resindex=atom_resindex,
            trajectory=True,
        )
        concatenated_traj.add_TopologyAttr("resname", [self.lipid_resname])
        concatenated_traj.load_new(coords, in_memory=False)

        return concatenated_traj, n_lipid, n_frames * n_lipid

    def concatenate_traj_with_merging(self):
        """
        ConcatenateTrajWithMerging performs basic trajectory concatination. In
        contrast to basic concatenateTraj it additionally merges splitted lipids.
        First, it creates extracts coordinates from trajectory, next, it reshapes
        the coordinate array, swaps time and resid axes to obtain continuous
        trajectories of individual lipids (this is needed for autocorrelation
        time analysis), and finally merges individual lipid trajectories.
        """
        traj = self.traj.trajectory
        n_frames = len(traj)

        heavy_atoms_topology = self.headlist[0] + self.tail1list[0] + self.tail2list[0]

        for head, sn1, sn2 in zip(
            self.headlist[1:],
            self.tail1list[1:],
            self.tail2list[1:],
            strict=False,
        ):
            heavy_atoms_topology = heavy_atoms_topology.union(head + sn1 + sn2)

        n_atoms_lipid = len(self.topology.atom_names())

        n_lipid = heavy_atoms_topology.n_atoms // n_atoms_lipid

        # TODO: add check
        # n_atoms_lipid == heavy_atoms_topology.n_atoms

        # Get all coordinates n_frames, n_lipid * n_atoms_lipid
        coords = AnalysisFromFunction(lambda ag: ag.positions.copy(), heavy_atoms_topology).run().results["timeseries"]
        # Adding axis to separate individual lipid trajectories
        coords = coords.reshape((n_frames, n_lipid, n_atoms_lipid, 3))
        # Swapping time frame with lipid axis
        # Now we have continuous lipid trajectory
        coords = np.swapaxes(coords, 0, 1)
        # Reshaping to desired shape
        coords = coords.reshape((n_frames * n_lipid, n_atoms_lipid, 3))

        # Creating new Universe from the concatenated coordinates
        atom_resindex = np.repeat(0, n_atoms_lipid)
        concatenated_traj = mda.Universe.empty(
            n_atoms_lipid,
            1,
            atom_resindex=atom_resindex,
            trajectory=True,
        )
        concatenated_traj.add_TopologyAttr("resname", [self.lipid_resname])
        concatenated_traj.load_new(coords, in_memory=False)

        return concatenated_traj, n_lipid, n_frames * n_lipid

    def align_traj(self, concatenated_traj):
        """
        AlignTraj alignes the concatenated trajectory in two steps: (1) it computes
        average structure after alignment to the first frame, and (2) it alignes
        the structure to the calculated average structure in (1).
        """
        # Compute average structure after alignment to the first frame
        av = align.AverageStructure(concatenated_traj, ref_frame=0).run()
        # Align to average structure
        align.AlignTraj(concatenated_traj, av.results.universe).run()
        # Compute average structure after second alignment
        coords = (
            AnalysisFromFunction(
                lambda ag: ag.positions.copy(),
                concatenated_traj.select_atoms("all"),
            )
            .run()
            .results["timeseries"]
        )
        return coords.reshape(-1, 3), coords.mean(axis=0).reshape(1, -1)

    def concatenate(self):
        """Simple enveloping function to perform concatenation."""
        print(f"Concatenator: Concatenating lipid with resname {self.lipid_resname}")
        if not self.topology.is_merge_needed():
            # Merging is not needed
            concatenated_traj, n_lipid, n_frames = self.concatenate_traj()
        else:
            concatenated_traj, n_lipid, n_frames = self.concatenate_traj_with_merging()
        aligned_traj, av_pos = self.align_traj(concatenated_traj)

        return aligned_traj, av_pos, n_lipid, n_frames


class PCA:
    """
    Class PCA is a class that actually performs PCA. It has the following methods:

    1. Simple constructor, which sets the aligned trajtory, average coordinates,
       number of lipids, number of frames in the concatenated trajectory and
       trajectory length in ns
    2. :meth:`PCA` prepares the trajectory coordinates for analysis and calculates
       principal components
    3. :meth:`get_proj` projects the trajectory on principal components
    4. :meth:`get_lipid_autocorrelation` calculates the autocorrelation timeseries for
       individual lipid
    5. :meth:`get_autocorrelations` calculates the autocorrelation timeseries for
       trajectory

    :param aligned_traj: np.array with positions of concatenated and aligned
                         trajectory
    :param av_pos: np.array with average positions for the lipid
    :param n_lipid: number of lipids of particular type in the system
    :param n_frames: number of frames in the concatenated trajectory
    :param traj_time: trajectory length in ns
    """

    def __init__(self, aligned_traj, av_pos, n_lipid, n_frames, traj_time):
        self.aligned_traj = aligned_traj
        self.av_pos = av_pos
        self.n_lipid = n_lipid
        self.n_frames = n_frames
        self.traj_time = traj_time

    def PCA(self):  # noqa: N802
        """
        PCA calculates the PCA. First the data is centered and then covariance
        matrix is calculated manually.
        """
        # centering of positions relative to the origin
        x = self.aligned_traj.astype(np.float64)
        x = x.reshape(self.n_frames, self.av_pos.shape[1]) - self.av_pos
        # the sum of all coordinates (to calculate mean)
        x_1 = x.sum(axis=0)
        # production of X and X^T
        x_x = np.tensordot(x, x, axes=(0, 0))
        # covariance matrix calculation
        cov_mat = (x_x - np.dot(x_1.reshape(len(x_1), 1), (x_1.reshape(len(x_1), 1)).T) / self.n_frames) / (
            self.n_frames - 1
        )
        # eigenvalues and eigenvectors calculation
        eig_vals, eig_vecs = np.linalg.eigh(cov_mat)
        self.eig_vecs = np.flip(eig_vecs, axis=1).T

        return x

    def get_proj(self, cdata):
        """Projecting the trajectory on the 1st principal component."""
        # projection on PC1
        proj = np.tensordot(cdata, self.eig_vecs[0:1], axes=(1, 1)).T
        proj = np.concatenate(np.array(proj), axis=None)

        self.proj = proj

    def get_lipid_autocorrelation(self, data, variance, mean):
        """Autocorrelation calculation for individual lipid."""
        # Centering the data
        data -= mean
        # Convolve the data
        r = signal.fftconvolve(data, data[::-1], mode="full")[-len(data) :]
        # Weight the correlation to get the result in range -1 to 1
        return r / (variance * (np.arange(len(data), 0, -1)))

    def get_autocorrelations(self):
        """Autocorrelation calculation for the trajectory."""
        variance = self.proj.var()
        mean = self.proj.mean()
        # extract the trajectories for individual lipids
        separate_projs = [
            self.proj[len(self.proj) * i // self.n_lipid : len(self.proj) * (i + 1) // self.n_lipid]
            for i in range(self.n_lipid)
        ]
        # calculate autocorrelations for individual lipids
        r = np.array(
            [self.get_lipid_autocorrelation(x, variance, mean) for x in separate_projs],
        )
        r = r.mean(axis=0)
        t = np.arange(len(r)) * self.traj_time / len(r)
        self.autocorrelation = np.array([t, r]).T


class TimeEstimator:
    """
    Class TimeEstimator is a class that estimates equilbration time from
    autocorrelation data.
    It includes the following methods:

    1. Simple constructor, which sets the autocorrelation data
    2. Helper method :meth:`get_nearest_value` that finds the interval in the data
       where the autocorrelation falls below specific value
    3. :meth:`timerelax` method calculates the logarithms of autocorrelation and
       calculates the decay time
    4. :meth:`calculate_time` method calculates the estimated equilibration time

    :param autocorrelation: autocorrelation data
    """

    def __init__(self, autocorrelation):
        self.autocorrelation = autocorrelation

    def get_nearest_value(self, iterable, value):
        """
        get_nearest_value method return the indices that frame the particular value
        As an input it gets an array, which is expected to decay. The method tries
        to find and index (index_A), for which the array gets below the cutoff value
        for the first time. Then, for index_A-1 the array was larger than the cutoff
        value for the last time. It means that the function, represented by the
        array data is equal to cutoff somewhere between timesteps that correspond to
        index_A-1 and index_A. It can happen, that this index_A is equal 0. Then,
        method searches for the first occurance, where array is smaller than
        array[0]. The method returns the interval inbetween the array gets below
        cutoff.
        """
        try:
            a = np.where(iterable < value)[0][0]
        except Exception:
            print("TimeEstimator: Autocorrelations do not converge. We shift to extrapolation regime.")
            a = np.where(iterable == np.min(iterable))[0][0]
            return a, a - 1

        if a == 0:
            b = np.where(iterable < iterable[a])[0][0]
        else:
            b = a - 1
        return a, b

    def timerelax(self):
        """Estimates autocorrelation decay time."""
        time = self.autocorrelation[:, 0]
        autocorrelation = self.autocorrelation[:, 1]

        points = []
        j = 1
        while j < len(autocorrelation):
            points.append(j)
            j = int(1.5 * j) + 1
        # Get timeseries for autocorrelation
        t_pic = np.array([time[i] for i in points if autocorrelation[i] > 0.0])
        r_pic = np.array(
            [autocorrelation[i] for i in points if autocorrelation[i] > 0.0],
        )
        # Calculate logs for time and autocorrelation arrays
        # We use log since in log-log scale autocorrelation is closer to linear
        r_log = np.log(r_pic)
        t_log = np.log(t_pic)
        # data interpolation. We are searching for time interval where
        # autocorrelations decay by e. This is the most stable method.
        # Autocorrelations decay by e is equivalent to
        # log(autocorrelation) < - 1
        power = -1
        ai, bi = self.get_nearest_value(r_log, power)
        # perform interpolation
        a = (t_log[bi] - t_log[ai]) / (r_log[bi] - r_log[ai])
        b = t_log[ai] - a * r_log[ai]
        t_relax1 = a * power + b
        t_relax1 = np.exp(t_relax1)

        return t_relax1

    def calculate_time(self):
        """
        Basic enveloping method that estimates equilibration time from
        autocorrelation decay time. They are linearly connected and the
        coefficient is calculated experimentally.
        """
        # relaxation time at e^1 decay
        te1 = self.timerelax()

        return te1 * 49  # 49 from sample data average of tkss2/tac1
