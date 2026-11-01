"""Tests for loading experiments used by experiment matching."""

import pytest

pytestmark = [pytest.mark.min, pytest.mark.sim2]


def test_experiments_are_readable_through_package_api():
    """ToyData experiments are discoverable through the public collection API."""
    from fairmd.lipids.experiment import ExperimentCollection

    ExperimentCollection.clear_instance()
    op_experiments = ExperimentCollection.load_from_data("OPExperiment")
    assert len(op_experiments) == 1
    assert op_experiments[0].exp_id == "toy_op"
    assert op_experiments[0].data["POPC"]["M_G1_M M_G1H1_M"][:2] == [-0.17, 0.02]

    ExperimentCollection.clear_instance()
    ff_experiments = ExperimentCollection.load_from_data("FFExperiment")
    assert len(ff_experiments) == 1
    assert ff_experiments[0].exp_id == "toy_ff"
    assert ff_experiments[0].data[1] == [0.05, 0.80, 0.04]


def test_match_experiments_matching_criteria():
    """Matching criteria accept compatible simulation and experiment values."""
    from fairmd.lipids.matcher import E2SComparator

    assert E2SComparator.match_membrane_composition({"POPC": 1.0}, {"POPC": 1.0})
    assert E2SComparator.match_solution_composition({"SOD": 0.10}, {"SOD": 0.12})
    assert E2SComparator.match_temperature(303.0, 303.15)
    assert not E2SComparator.match_temperature(303.0, 310.0)


def test_find_pairs_and_change_sims_works(tmp_path):
    """Matching pairs are returned and written to the simulation metadata."""
    from fairmd.lipids.matcher import pair_all_simulations

    class FakeSimulation(dict):
        def __init__(self, simulation_id, path, temperature=303.0, implicit_water=False):
            super().__init__(
                ID=simulation_id,
                TEMPERATURE=temperature,
                EXPERIMENT={"ORDERPARAMETER": {}, "FORMFACTOR": []},
                path=str(path),
            )
            self.implicit_water = implicit_water

        @property
        def readme(self):
            return {key: value for key, value in self.items()}

        def to_readme(self):
            return self.readme

        def membrane_composition(self, basis):
            assert basis == "molar"
            return {"POPC": 1.0}

        def solution_composition(self, basis):
            assert basis == "molar"
            if self.implicit_water:
                raise ValueError("implicit water")
            return {"SOD": 0.1}

        def get_hydration(self, basis):
            assert basis == "number"
            return 20.0

    class FakeExperiment(dict):
        def __init__(self, exp_id, exptype, data, temperature=303.0):
            super().__init__(TEMPERATURE=temperature)
            self.exp_id = exp_id
            self.exptype = exptype
            self.data = data

        def membrane_composition(self, basis):
            return {"POPC": 1.0}

        def solution_composition(self, basis):
            return {"SOD": 0.1}

        def get_hydration(self, basis):
            return 20.0

    simulation = FakeSimulation("matching", path=tmp_path)
    skipped = FakeSimulation("implicit", path=tmp_path, implicit_water=True)
    op1_experiment = FakeExperiment("op-1", "OrderParameters", {"POPC": []}, temperature=303)
    op2_experiment = FakeExperiment("op-2", "OrderParameters", {"POPC": []}, temperature=270)

    with pytest.warns(UserWarning, match="implicit water"):
        pairs = pair_all_simulations([op1_experiment, op2_experiment], [simulation, skipped])

    assert pairs == [[simulation, op1_experiment]]
    assert simulation["EXPERIMENT"] == {
        "ORDERPARAMETER": {"POPC": ["op-1"]},
        "FORMFACTOR": [],
    }
    assert skipped["EXPERIMENT"] == {"ORDERPARAMETER": {}, "FORMFACTOR": []}


def test_match_experiments_integration(tmp_path, monkeypatch):
    """The matcher loads ToyData and writes its pair-search log."""
    from fairmd.lipids.bin.match_experiments import match_experiments

    monkeypatch.chdir(tmp_path)
    match_experiments()

    log = (tmp_path / "search-databank-pairs.log").read_text()
    assert "=== OP PAIRS ===" in log
    assert "=== FF PAIRS ===" in log

    from fairmd.lipids.core import initialize_databank

    ss = initialize_databank()
    s = ss.loc(281)

    assert s["EXPERIMENT"]["ORDERPARAMETER"] == {"POPC": ["toy_op"]}
    assert s["EXPERIMENT"]["FORMFACTOR"] == ["toy_ff"]
