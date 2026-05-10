"""
Accumulate functions required for major evaluate_quality script.

TODO: add proper tests
"""

import json
import os
import warnings
from abc import ABC, abstractmethod

import numpy as np
import numpy.typing as npt
import scipy.stats

import fairmd.lipids.analib.formfactor as ff
from fairmd.lipids import FMDL_SIMU_PATH
from fairmd.lipids.api import get_FF, get_OP
from fairmd.lipids.auxiliary import CompactJSONEncoder, mollib
from fairmd.lipids.core import System, initialize_databank
from fairmd.lipids.experiment import ExperimentCollection, OPExperiment


class QualSimulation(System):
    def __init__(self, s: System):
        super().__init__(s)
        self.op_data = get_OP(s)
        try:
            self.ff_data = get_FF(s)
        except FileNotFoundError:
            msg = f"System #{s['ID']} doesn't have computed FF data."
            warnings.warn(msg, stacklevel=2)
            self.ff_data = None

    @staticmethod
    def load_all_paired(id_list: list[int] | None = None) -> list["QualSimulation"]:
        """Load simulations with experimental pairings."""
        systems = initialize_databank()

        simulations = []
        for system in systems:
            if id_list is not None and system["ID"] not in id_list:
                continue
            experiments = system.get("EXPERIMENT", {})
            if any(experiments.values()):  # if experiments is not empty
                s = QualSimulation(system)
                simulations.append(s)
        print(f"Loaded {len(simulations)} that has some experimental pairings.")

        return simulations


class QualityEvaluator(ABC):
    """Interface for quality evaluation of simulation."""

    def __init__(self, simulation: QualSimulation, experiments: ExperimentCollection):
        self._sim = simulation
        self._exps = experiments

    @abstractmethod
    def evaluate_one(self) -> bool:
        """Evaluate quality for one simulation."""

    @abstractmethod
    def save_results(self) -> None:
        """Save quality evaluation results."""

    @staticmethod
    def prob_2_within_trustinterval(
        xv: npt.ArrayLike,
        xerr: npt.ArrayLike,
        yv: npt.ArrayLike,
        yerr: npt.ArrayLike,
    ) -> npt.ArrayLike:
        """
        Compute the probability of two vals being a part of one distribution.

        Used quality value from experimental and simulation OP data.
        Probability is computed using Eq. (3) [10.1038/s41467-024-45189-z].

        NOTE: Computing the probability taking into account small values using
        scipy.special.log1p is not required if sd is above 1e-5, which is
        currently the case for all OP data. It is also not required if we are
        not going to bring values to log-scale.

        Args:
            xv: Experimental OP value (all float/arrays)
            xerr: Experimental error
            yv: Simulated OP value
            yerr: Standard deviation from simulation

        Returns
        -------
            float/ndarray: probability value(s) or nans
        """
        # normal distribution N(s, OP_sim, op_sim_sd)
        a = xv - xerr
        b = xv + xerr

        a_rel = (yv - a) / yerr
        b_rel = (yv - b) / yerr

        p_b = scipy.stats.t.sf(b_rel, df=1, loc=0, scale=1)
        p_a = scipy.stats.t.sf(a_rel, df=1, loc=0, scale=1)

        return p_b - p_a

    @staticmethod
    def _round_quality_values(obj: dict | list, ndigits: int = 4) -> dict | list:
        """Round all floating-point values in a nested dict/list structure."""
        stack = [obj]

        while stack:
            current = stack.pop()

            if isinstance(current, dict):
                for k, v in current.items():
                    if isinstance(v, float):
                        current[k] = round(v, ndigits)
                    elif isinstance(v, (dict, list)):
                        stack.append(v)

            elif isinstance(current, list):
                for i, v in enumerate(current):
                    if isinstance(v, float):
                        current[i] = round(v, ndigits)
                    elif isinstance(v, (dict, list)):
                        stack.append(v)

        return obj


class OPQualityEvaluator(QualityEvaluator):
    """Evaluate quality of order parameter data."""

    def evaluate_one(self) -> bool:
        self.lipid_fragment_qualities = {}
        self.lipid_atomic_qualities = {}

        for lipidname, lipid in self._sim.lipids.items():
            md_lipid_ops = self._sim.op_data[lipidname]
            if md_lipid_ops is None:
                continue

            # TODO: bb is merged into headgroup. But sn-s do not.. Cryptic rule.
            # TODO: What to do with other types of lipids?
            fragments = mollib.get_fragments(lipid.mapping_dict)
            fragment_qual_perexp = {}
            lipid_quality_perexp = {}

            for expid in self._sim["EXPERIMENT"]["ORDERPARAMETER"].get(lipidname, []):
                print(f"OP quality of simulation data in {self._sim['path']}")
                print(
                    f".. evaluating {lipidname} lipid using experimental data from {expid}",
                )
                exp_lipid_ops = self._exps.loc(expid).data[lipidname]
                # save qualities of simulation-vs-experiment into a dictionary
                lipid_quality_perexp[expid] = self.atomic_quality(exp_lipid_ops, md_lipid_ops)

                # calculate quality for molecule fragments headgroup, sn-1, sn-2
                _frq = self.atomic2fragment_quality(lipid_quality_perexp[expid], fragments)
                _frw = self.get_fragments_coverage(fragments, lipid_quality_perexp[expid])
                # dot product [ qualities * weights ]
                fragment_qual_perexp[expid] = {k: _frq[k] * _frw[k] for k in fragments}

            if not fragment_qual_perexp:
                continue

            fragment_quality_merged = self.fragment_quality_unite_multexp(fragment_qual_perexp)

            # export per-lipid qualities
            self.lipid_fragment_qualities[lipidname] = fragment_quality_merged
            self.lipid_atomic_qualities[lipidname] = lipid_quality_perexp

        if not self.lipid_fragment_qualities:
            return False

        self.system_qual_output = self.system_quality_gather_lipids()
        return True

    def save_results(self) -> None:
        wdir = os.path.join(FMDL_SIMU_PATH, self._sim["path"])

        # write into the XXXOrderParameters_quality.json quality data file for each lipid
        for lipidname, atomic_quality_perexp in self.lipid_atomic_qualities.items():
            fname = os.path.join(wdir, lipidname + "_OrderParameters_quality.json")
            self._round_quality_values(atomic_quality_perexp)
            with open(fname, "w") as f:
                json.dump(atomic_quality_perexp, f, cls=CompactJSONEncoder)

        # Write fragment qualities for each lipid
        for lipidname, fragment_quality_merged in self.lipid_fragment_qualities.items():
            if any(v > 0 for v in fragment_quality_merged.values()):
                fname = os.path.join(wdir, lipidname + "_FragmentQuality.json")
                self._round_quality_values(fragment_quality_merged)
                with open(fname, "w") as f:
                    json.dump(fragment_quality_merged, f)

        # Write system quality
        if any(v > 0 for v in self.system_qual_output.values()):
            fname = os.path.join(wdir, "SYSTEM_quality.json")
            self._round_quality_values(self.system_qual_output)
            with open(fname, "w") as f:
                json.dump(self.system_qual_output, f)

        print("Order parameter quality written " + self._sim["path"])

    @classmethod
    def atomic_quality(cls, exp_op_data: dict, sim_op_data: dict) -> dict:
        """
        Calculate quality for a molecule (times their weights in exp data).

        :param exp_op_data: dictionary of type {op_uname: [op_value]}.
        :param sim_op_data: dictionary of type {op_uname: [op_value, op_sigma, op_sd]}.

        :return: dictionary of type {"nC nH": quality value}.
        """
        # TODO: hardcoded error value, should be taken from experiment data when available
        exp_error = OPExperiment.DEFAULT_ERROR

        # union of keys in exp_op_data and sim_op_data
        all_keys = sorted(set(exp_op_data.keys()) | set(sim_op_data.keys()))
        res_dict = {}

        for key in all_keys:
            if key not in exp_op_data or key not in sim_op_data:
                res_dict[key] = np.nan
                continue
            q = cls.prob_2_within_trustinterval(
                xv=exp_op_data[key][0],
                xerr=exp_error,
                yv=sim_op_data[key][0],
                yerr=sim_op_data[key][2],
            )
            res_dict[key] = q

        return res_dict

    @classmethod
    def get_fragments_coverage(cls, fragments: dict, q_data: dict) -> dict:
        """
        Calculate the coverage in data of each fragment.

        The coverage is implemented as the percent of non-nans.

        :param fragments: Dictionary of type {fragment: lists of unames}.
        :param q_data: Dictionary of type {op_uname: quality value}.

        :return: Dictionary of type {fragment: percentage of evaluated OPs}
        """
        frag_percentage = dict.fromkeys(fragments, 0)

        for frg_name, frg_atoms in fragments.items():
            frg_atoms_set = set(frg_atoms)
            nan_mask = [np.isnan(v) for k, v in q_data.items() if k.split(" ")[0] in frg_atoms_set]
            frag_percentage[frg_name] = 1 - sum(nan_mask) / len(nan_mask) if len(nan_mask) > 0 else 0

        return frag_percentage

    @classmethod
    def atomic2fragment_quality(cls, atomic_qual_dict: dict, fragments: dict) -> dict:
        fqdict = dict.fromkeys(fragments.keys())
        for frg_name, frg_atoms in fragments.items():
            q_list = [v for k, v in atomic_qual_dict.items() if k.split(" ")[0] in frg_atoms]
            fqdict[frg_name] = np.nanmean(q_list) if len(q_list) > 0 else np.nan
        return fqdict

    @classmethod
    def fragment_quality_unite_multexp(cls, fragment_qual_dict: dict) -> dict:
        """
        Condition fragment qualities.

        The second-layer function.

        :param fragment_qual_perexp: dictionary of type {expid: {fragment: quality value}}.

        :return: dictionary of type {fragment: average quality value}.
        """
        sums_dict = {}
        for exp_frag_qual in fragment_qual_dict.values():
            for frag_name, frag_qual in exp_frag_qual.items():
                sums_dict.setdefault(frag_name, []).append(frag_qual)

        avg_total_quality = {frag_name: np.nanmean(frag_vals) for frag_name, frag_vals in sums_dict.items()}
        avg_total_quality["total"] = np.nanmean(list(avg_total_quality.values()))

        return avg_total_quality

    def system_quality_gather_lipids(self) -> dict:
        """
        Gather fragment qualities for each lipid in a system and compute system-wide quality.

        :return: dictionary of type {macrofragment: quality} where macrofragmetn is headgroup|tails|total.
        """
        system_fragment_qualities = self.lipid_fragment_qualities
        molar_composition = self._sim.membrane_composition(basis="molar")

        system_dict = {}
        lipid_dict = {}
        w_nan = []

        for lname, lqual in system_fragment_qualities.items():
            lipid_dict = dict.fromkeys(lqual.keys(), 0)
            w = molar_composition[lname]
            for key, value in lqual.items():
                if not np.isnan(value):
                    lipid_dict[key] += w * value
                else:
                    # save 1 - w of a lipid into a list if the fragment quality is nan
                    w_nan.append(1 - w)
            system_dict[lname] = lipid_dict

        system_quality = {}

        headgroup = 0
        tails = 0
        total = 0

        for ldict in system_dict.values():
            for key, value in ldict.items():
                if key == "total":
                    total += value
                elif key == "headgroup":
                    headgroup += value
                elif key in {"sn-1", "sn-2"}:
                    tails += value / 2
                else:
                    tails += value  # everything non head is tail??

        penalty = np.prod(w_nan) if w_nan else 1

        system_quality["headgroup"] = headgroup * penalty
        system_quality["tails"] = tails * penalty
        system_quality["total"] = total * penalty

        return system_quality


class FFQualityEvaluator(QualityEvaluator):
    """Evaluate quality of form factor data."""

    def evaluate_one(self, method: str = "new") -> bool:
        results_ff = {}
        if method == "new":
            method_func = self.calc_new_ff_quality
        elif method == "old":
            method_func = self.calc_ff_quality

        if "FORMFACTOR" not in self._sim.get("EXPERIMENT", {}) or not self._sim["EXPERIMENT"]["FORMFACTOR"]:
            return False
        if self._sim.ff_data is None:
            return False
        for expid in self._sim["EXPERIMENT"]["FORMFACTOR"]:
            exp_ff_data = np.array(self._exps.loc(expid).data)
            results_ff[expid] = [
                method_func(self._sim.ff_data, exp_ff_data),
                ff.calc_ff_scaling_distance(exp_ff_data, self._sim.ff_data)[0],
            ]

        # TODO: handle multiple FF experiments better
        # currently, just pick the best one
        best_ep = None
        for expid, ffq_scf in results_ff.items():
            if ffq_scf is None:
                return False
            if best_ep is None or ffq_scf[0] < results_ff[best_ep][0]:
                best_ep = expid

        print(f"Form factor quality used for experiment data from {best_ep}")
        self.ff_quality_output = results_ff[best_ep]

        return True

    def save_results(self) -> None:
        print("Distance =", self.ff_quality_output[0], "; scaling factor =", self.ff_quality_output[1])

        wdir = os.path.join(FMDL_SIMU_PATH, self._sim["path"])
        print("Form factor quality saved for ", wdir)
        fname = os.path.join(wdir, "FormFactorQuality.json")

        ff_quality = list(self.ff_quality_output)
        self._round_quality_values(ff_quality)

        with open(fname, "w") as f:
            json.dump(ff_quality, f)

    @staticmethod
    def calc_ff_quality(ffd_sim: np.ndarray, ffd_exp: np.ndarray) -> float:
        """
        Calculate form factor quality.

        Quality calculation is performed as defined by Kučerka et al. 2010, doi:10.1007/s00232-010-9254-5

        :param ffd_sim: Simulation FF data (float 2D list)
        :param ffd_exp: Experiment FF data (float 2D list)

        :return: Quality value.
        """
        sim_min = ff.get_mins_from_ffdata(ffd_sim)
        exp_min = ff.get_mins_from_ffdata(ffd_exp)

        return np.abs(sim_min[0] - exp_min[0]) * 100

    @classmethod
    def calc_new_ff_quality(cls, ffd_sim: np.ndarray, ffd_exp: np.ndarray) -> float:
        """Calculate form factor quality via p-value of the difference between minima."""
        sim_pos, sim_err = ff.calc_minpos_with_error(ffd_sim, 0.5)
        exp_pos, exp_err = ff.calc_minpos_with_error(ffd_exp, 0.05)
        print(f"Sim: {sim_pos}+-{sim_err}; Exp: {exp_pos}+-{exp_err}")

        z = abs(sim_pos - exp_pos) / np.sqrt(sim_err**2 + exp_err**2)
        dpos_ref = 0.02  # natural physical scale
        sigma_ref = 0.01
        consistency_weight = 0.3

        # Term 1: are the positions close? (penalises large absolute difference)
        distance_score = np.exp(-((sim_pos - exp_pos) ** 2) / dpos_ref**2)

        # Term 2: are the measurements precise? (penalises large uncertainty)
        precision_score = np.exp(-(sim_err**2 + exp_err**2) / 4 * sigma_ref**2)

        # Term 3: are they statistically consistent? (penalises large z)
        consistency_score = np.exp(-consistency_weight * z**2)

        return distance_score * precision_score * consistency_score
