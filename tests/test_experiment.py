"""
:module: test_experiment.py

:description: Unittest for the experiment module
"""

import json
import os
import numpy as np
import pytest
import pytest_check as check
import yaml


pytestmark = [pytest.mark.min, pytest.mark.nodata]


@pytest.fixture
def mock_op_experiment_path(tmpdir):
    """Fixture for a mock order parameter experiment path."""
    exp_dir = tmpdir.mkdir("OPExperiment_1")
    readme_content = {
        "MOLAR_FRACTIONS": {"POPC": 1.0},
        "TEMPERATURE": 303.15,
        "ION_CONCENTRATIONS": {"NA": 0.1},
        "COUNTER_IONS": ["CL"],
    }
    with open(exp_dir.join("README.yaml"), "w") as f:
        yaml.dump(readme_content, f)

    data_content = {"C22": [-0.1, -0.2], "C23": [-0.15, -0.25]}
    with open(exp_dir.join("POPC_OrderParameters.json"), "w") as f:
        json.dump(data_content, f)

    return str(exp_dir)


@pytest.fixture
def mock_ff_experiment_path(tmpdir):
    """Fixture for a mock form factor experiment path."""
    exp_dir = tmpdir.mkdir("FFExperiment_1")
    readme_content = {
        "MOLAR_FRACTIONS": {"DOPC": 1.0},
        "TEMPERATURE": 298.15,
    }
    with open(exp_dir.join("README.yaml"), "w") as f:
        yaml.dump(readme_content, f)

    data_content = {"q": [0.1, 0.2, 0.3], "I": [1.0, 0.5, 0.2]}
    with open(exp_dir.join("system_FormFactor.json"), "w") as f:
        json.dump(data_content, f)

    return str(exp_dir)


@pytest.fixture
def mock_missing_readme_path(tmpdir):
    """Fixture for a mock experiment path missing README.yaml."""
    return str(tmpdir.mkdir("Missing_Readme"))


@pytest.fixture
def mock_empty_data_path(tmpdir):
    """Fixture for a mock experiment path with no data files."""
    exp_dir = tmpdir.mkdir("Empty_Data")
    readme_content = {"TEMPERATURE": 300.0}
    with open(exp_dir.join("README.yaml"), "w") as f:
        yaml.dump(readme_content, f)
    return str(exp_dir)


class TestOPExperiment:
    """Test the OPExperiment class."""

    def test_initialization(self, mock_op_experiment_path):
        from fairmd.lipids.experiment import OPExperiment

        """Test successful initialization."""
        exp = OPExperiment("exp1", mock_op_experiment_path)
        assert exp.exp_id == "exp1"
        assert exp.path == mock_op_experiment_path
        assert exp.metadata["TEMPERATURE"] == 303.15

    def test_missing_readme_raises_error(self, mock_missing_readme_path):
        from fairmd.lipids.experiment import OPExperiment

        """Test that FileNotFoundError is raised if README.yaml is missing."""
        with pytest.raises(FileNotFoundError):
            OPExperiment("exp_fail", mock_missing_readme_path)

    def test_data_loading(self, mock_op_experiment_path):
        from fairmd.lipids.experiment import OPExperiment

        """Test correct loading of order parameter data."""
        exp = OPExperiment("exp1", mock_op_experiment_path)
        data = exp.data
        assert "POPC" in data
        assert data["POPC"]["C22"] == [-0.1, -0.2]

    def test_properties(self, mock_op_experiment_path):
        from fairmd.lipids.experiment import OPExperiment

        """Test properties of OPExperiment."""
        exp = OPExperiment("exp1", mock_op_experiment_path)
        check.equal(exp.exptype, "OrderParameters")
        check.equal(OPExperiment.target_folder(), "OrderParameters")

    def test_empty_data(self, mock_empty_data_path):
        from fairmd.lipids.experiment import OPExperiment

        """Test behavior with no data files."""
        exp = OPExperiment("exp_empty", mock_empty_data_path)
        assert exp.data == {}


class TestFFExperiment:
    """Test the FFExperiment class."""

    def test_initialization(self, mock_ff_experiment_path):
        from fairmd.lipids.experiment import FFExperiment

        """Test successful initialization."""
        exp = FFExperiment("exp2", mock_ff_experiment_path)
        assert exp.exp_id == "exp2"
        assert exp.path == mock_ff_experiment_path
        check.almost_equal(exp.metadata["TEMPERATURE"], 298.15, abs=1e-6)

    def test_data_loading(self, mock_ff_experiment_path):
        from fairmd.lipids.experiment import FFExperiment

        """Test correct loading of form factor data."""
        exp = FFExperiment("exp2", mock_ff_experiment_path)
        data = exp.data
        assert "q" in data
        np.testing.assert_allclose(data["I"], [1.0, 0.5, 0.2])

    def test_properties(self, mock_ff_experiment_path):
        from fairmd.lipids.experiment import FFExperiment

        """Test properties of FFExperiment."""
        exp = FFExperiment("exp2", mock_ff_experiment_path)
        assert exp.exptype == "FormFactors"
        assert FFExperiment.target_folder() == "FormFactors"

    def test_empty_data(self, mock_empty_data_path):
        from fairmd.lipids.experiment import FFExperiment

        """Test behavior with no data files."""
        exp = FFExperiment("exp_empty", mock_empty_data_path)
        assert exp.data == {}


class TestExperimentBase:
    """Test base Experiment class functionality through a subclass."""

    def test_get_lipids(self, mock_op_experiment_path):
        from fairmd.lipids.experiment import OPExperiment

        """Test the get_lipids method."""
        exp = OPExperiment("exp1", mock_op_experiment_path)
        lipids = exp.get_lipids()
        assert "POPC" in lipids
        assert "DOPC" not in lipids

    def test_get_ions(self, mock_op_experiment_path):
        from fairmd.lipids.experiment import OPExperiment

        """Test the get_ions method."""
        exp = OPExperiment("exp1", mock_op_experiment_path)
        ions = exp.get_ions(ions=["NA", "CL", "K"])
        assert "NA" in ions
        assert "CL" in ions
        assert "K" not in ions

    def test_dunder_methods(self, mock_op_experiment_path):
        from fairmd.lipids.experiment import OPExperiment

        """Test __repr__, __eq__, __hash__, and __getitem__."""
        exp1 = OPExperiment("exp1", mock_op_experiment_path)
        exp2 = OPExperiment("exp1", mock_op_experiment_path)
        exp3 = OPExperiment("exp3", mock_op_experiment_path)

        assert repr(exp1) == "OPExperiment(id='exp1')"
        assert exp1 == exp2
        assert exp1 != exp3
        assert hash(exp1) == hash(exp2)
        assert hash(exp1) != hash(exp3)
        assert exp1["TEMPERATURE"] == 303.15

    def test_abstract_methods_not_implemented(self):
        from fairmd.lipids.experiment import Experiment

        """Test that abstract methods raise NotImplementedError."""
        # We can't instantiate Experiment directly, so we check a class method
        with pytest.raises(NotImplementedError):
            Experiment.target_folder()


@pytest.fixture
def mock_exp_data_path(tmpdir, monkeypatch):
    """Fixture to create a mock data structure and monkeypatch FMDL_EXP_PATH."""
    mock_root = tmpdir.mkdir("fairmd_data")
    op_path = mock_root.mkdir("OrderParameters").mkdir("exp1")
    ff_path = mock_root.mkdir("FormFactors").mkdir("exp2")

    # Create dummy OP experiment
    with open(op_path.join("README.yaml"), "w") as f:
        yaml.dump({"TEMPERATURE": 310.0}, f)
    with open(op_path.join("DUMMY_Order_Parameters.json"), "w") as f:
        f.write("{}")

    # Create dummy FF experiment
    with open(ff_path.join("README.yaml"), "w") as f:
        yaml.dump({"TEMPERATURE": 290.0}, f)
    with open(ff_path.join("system_FormFactor.json"), "w") as f:
        f.write("{}")

    # Monkeypatch the path to point to our mock directory
    monkeypatch.setattr("fairmd.lipids.experiment.FMDL_EXP_PATH", str(mock_root))
    return str(mock_root)


class TestExperimentCollection:
    """Test the ExperimentCollection class."""

    def test_singleton(self):
        from fairmd.lipids.experiment import ExperimentCollection

        """Test that ExperimentCollection is a singleton."""
        c1 = ExperimentCollection()
        c2 = ExperimentCollection()
        assert c1 is c2

    def test_add_and_get(self, mock_op_experiment_path):
        from fairmd.lipids.experiment import OPExperiment, ExperimentCollection

        """Test adding and retrieving an experiment."""
        ExperimentCollection.clear_instance()
        collection = ExperimentCollection()
        exp = OPExperiment("exp1", mock_op_experiment_path)
        collection.add(exp)
        assert collection.get("exp1") == exp
        assert "exp1" in collection

    def test_load_from_data(self, mock_exp_data_path):
        from fairmd.lipids.experiment import OPExperiment, FFExperiment, ExperimentCollection

        """Test loading experiments from a mock data directory."""
        ExperimentCollection.clear_instance()
        collection = ExperimentCollection.load_from_data()

        # Should load OP experiments by default
        assert len(collection) == 1
        assert "exp1" in collection
        assert isinstance(collection.get("exp1"), OPExperiment)
        assert collection.get("exp1").path == os.path.join(mock_exp_data_path, "OrderParameters", "exp1")

        # Test loading FF experiments
        ExperimentCollection.clear_instance()
        collection = ExperimentCollection()
        # Manually load FF experiments by iterating through types
        ff_path = os.path.join(mock_exp_data_path, FFExperiment.target_folder())
        for exp_id in os.listdir(ff_path):
            exp = FFExperiment(exp_id, os.path.join(ff_path, exp_id))
            collection.add(exp)

        assert len(collection) == 1
        assert "exp2" in collection
        assert isinstance(collection.get("exp2"), FFExperiment)

    def test_load_all_types(self, mock_exp_data_path):
        from fairmd.lipids.experiment import OPExperiment, FFExperiment, ExperimentCollection

        """Test loading all experiment types into one collection."""
        ExperimentCollection.clear_instance()
        collection = ExperimentCollection.load_from_data()  # Loads OP

        # Manually load FF and add to the same collection
        ff_path = os.path.join(mock_exp_data_path, FFExperiment.target_folder())
        for exp_id in os.listdir(ff_path):
            exp_path = os.path.join(ff_path, exp_id)
            if os.path.isdir(exp_path):
                exp = FFExperiment(exp_id, exp_path)
                collection.add(exp)

        assert len(collection) == 2
        assert "exp1" in collection
        assert "exp2" in collection
        assert isinstance(collection.get("exp1"), OPExperiment)
        assert isinstance(collection.get("exp2"), FFExperiment)

    def test_load_invalid_type_raises_error(self):
        from fairmd.lipids.experiment import ExperimentCollection

        """Test that loading an invalid experiment type raises ValueError."""
        ExperimentCollection.clear_instance()
        with pytest.raises(ValueError):
            ExperimentCollection.load_from_data(exp_type="InvalidType")
