"""Tests for the developer dataset-generation workflow."""

import logging
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest
import pytest_check as check

pytestmark = [pytest.mark.develop, pytest.mark.nodata]

DEVELOPER_DIR = Path(__file__).resolve().parents[2] / "developer"
sys.path.insert(0, str(DEVELOPER_DIR))


def test_cli_arguments(tmp_path, monkeypatch):
    """Smoke test of OP dataset generator CLI behavior"""
    import gen_op_dataset as cli

    check.equal(cli.main([]), 1, "gen_op_dataset requires arguments")

    calls = []

    def fake_generator(_: logging.Logger, hdf_fname: str) -> tuple[int, int]:
        calls.append(hdf_fname)
        return 1, 0

    monkeypatch.setattr(cli, "gen_op_from_exps", fake_generator)
    monkeypatch.setattr(cli, "gen_op_from_sims", fake_generator)
    monkeypatch.setattr(cli, "gen_opq_ds", fake_generator)

    expected_fname = {"exps": "exp-op-dataset.h5", "sims": "sims-op-dataset.h5", "qual": "atomic-opq-dataset.h5"}
    for arg0 in "exps", "sims", "qual":
        args = cli.build_parser().parse_args(["--" + arg0, "--output-dir", str(tmp_path)])
        calls = []
        result = cli.run(args)

        assert arg0 in result, f"requested {arg0} should be in results map"
        check.equal(result[arg0], (1, 0), "expected number of success and failed records differs")
        assert len(calls) == 1, "only one call is expected"

        check.equal(calls[0].rsplit("/", 1)[-1], expected_fname[arg0], "unexpected file stored")


class FakeSample(dict):
    """Small sample double implementing the interfaces used by dataset storers."""

    def __init__(self):
        super().__init__(ID="SIM-1", TEMPERATURE=300, COMPOSITION={"LIPID": {"COUNT": 2}})
        self.lipids = {
            "LIPID": SimpleNamespace(
                mapping_dict={"C1": {"SMILEIDX": 1}, "C2": {"SMILEIDX": 2}},
                metadata={"bioschema_properties": {"inChIKey": "KEY", "smiles": "C"}},
            ),
        }
        self.readme = {}

    def membrane_composition(self, basis: str) -> dict:
        assert basis == "molar"
        return {"LIPID": 1.0}

    def solution_composition(self, basis: str) -> dict:
        assert basis == "molar"
        return {"water": 1.0}

    def get_hydration(self):
        return 10.0


def test_common_op_conversion_sorts_smile_indices_and_values():
    from datasets.pureop import ExpOPDataStorer

    experiment = SimpleNamespace(
        exp_id="exp/1",
        data={"LIPID": {"C2 H": [0.2, 0.02], "C1 H": [0.1, 0.01]}},
        lipids=FakeSample().lipids,
        metadata={},
    )
    storer = ExpOPDataStorer(experiment, "LIPID")

    storer.prepare_dataframe()

    assert storer._df.to_dict("records") == [
        {"id": 1.0, "val": 0.1, "err": 0.01},
        {"id": 2.0, "val": 0.2, "err": 0.02},
    ]


def test_missing_smile_indices_are_rejected():
    from datasets.pureop import ExpOPDataStorer, OPDataError

    molecule = SimpleNamespace(mapping_dict={}, metadata={})
    experiment = SimpleNamespace(
        exp_id="exp-1", data={"LIPID": {"C1 H": [0.1]}}, lipids={"LIPID": molecule}, metadata={}
    )

    with pytest.raises(OPDataError, match="SMILEIDX"):
        ExpOPDataStorer(experiment, "LIPID").prepare_dataframe()


def test_quality_conversion_averages_three_equivalent_hydrogens():
    from datasets.opqual import OPQDataStorer

    storer = OPQDataStorer(FakeSample(), "LIPID")
    converted = storer._cvt_op_df(
        {"C1 H1": [0.1, 0.01], "C1 H2": [0.2, 0.02], "C1 H3": [0.3, 0.03]},
        err_pos=1,
    )

    row = converted.iloc[0]
    assert row["c"] == 1
    assert row["h"] == 1
    assert row["val"] == pytest.approx(0.2)
    assert row["err"] == pytest.approx(0.02)


def test_quality_points_inherit_simulation_sign():
    from datasets.opqual import OPQDataStorer

    storer = OPQDataStorer(FakeSample(), "LIPID")
    storer._sim_op = storer._cvt_op_df({"C1 H": [-0.4, 0.01]}, err_pos=1)
    storer._exp_opdict_one = storer._cvt_op_df({"C1 H": [0.2, 0.02]}, err_pos=1)

    storer.compute_q_points()

    assert storer._qpoints.iloc[0]["val_e"] == pytest.approx(-0.2)


def test_experiment_storer_writes_expected_hdf5_tables(tmp_path):
    from datasets.pureop import ExpOPDataStorer

    experiment = FakeSample()
    experiment.exp_id = "exp-1"
    experiment.data = {"LIPID": {"C1 H": [0.1, 0.01]}}
    experiment.metadata = {"NMR": {"METHOD": "ssNMR"}}
    experiment["TEMPERATURE"] = 298
    storer = ExpOPDataStorer(experiment, "LIPID")
    storer.prepare_dataframe()
    output = tmp_path / "dataset.h5"

    storer.store_to_hdf5(str(output))

    import pandas as pd

    with pd.HDFStore(output) as store:
        assert "/Eexp_1__LIPID/op_values" in store.keys()
        assert store["/Eexp_1__LIPID/op_values"].iloc[0]["val"] == pytest.approx(0.1)
