"""
Test quality module functions.

-------------------------------------------------------------------------------
NOTE: globally import of fairmd-lipids is **STRICTLY FORBIDDEN** because it
      breaks the substitution of global path folders
"""

import numpy as np
import numpy.testing as npt
import pytest
import pytest_check as check

# run only on sim2 mocking data
pytestmark = [pytest.mark.nodata, pytest.mark.min]


def test_prob_op_within_trustinterval():
    from fairmd.lipids.quality import QualityEvaluator

    # single value test
    op_exp = -0.22
    exp_error = 0.02
    op_sim = -0.21
    op_sim_sd = 0.017
    p1 = QualityEvaluator.prob_2_within_trustinterval(op_exp, exp_error, op_sim, op_sim_sd)
    check.almost_equal(0.5051486814438984, p1, abs=1e-6)

    # array test
    op_exp = np.array([-0.22, -0.18, -0.25])
    exp_error = np.array([0.02, 0.01, 0.03])
    op_sim = np.array([-0.21, -0.19, -0.23])
    op_sim_sd = np.array([0.017, 0.01, 0.02])

    p2 = QualityEvaluator.prob_2_within_trustinterval(op_exp, exp_error, op_sim, op_sim_sd)
    npt.assert_allclose(p2, [0.505148, 0.352416, 0.526465], atol=1e-6)

    # array test with nans
    op_exp = np.array([-0.22, -0.18, -0.25])
    exp_error = np.array([0.02, 0.01, 0.03])
    op_sim = np.array([-0.21, -0.19, -0.23])
    op_sim_sd = np.array([0.017, 0.0, 0.02])

    p2 = QualityEvaluator.prob_2_within_trustinterval(op_exp, exp_error, op_sim, op_sim_sd)
    npt.assert_allclose(p2, [0.505148, np.nan, 0.526465], atol=1e-6)


def test_calc_ff_quality_sin_exp_curve():
    from fairmd.lipids.quality import FFQualityEvaluator

    rng = np.random.default_rng()

    # Create FF-like data: |sin(x)| * exp(-x) over q-range typical for FF (0.1 to 3.0)
    q = np.linspace(0.005, 3.0, 1000)
    # Simulation data: minimum near q=0.5
    ffd_sim = np.column_stack([q, np.abs(np.sin(2 * np.pi * q)) * np.exp(-q)]).astype(float)

    # Test edge case: identical data
    quality_identical = FFQualityEvaluator.calc_ff_quality(ffd_sim, ffd_sim)
    check.almost_equal(0.0, quality_identical, abs=1e-6, msg="Identical data should yield zero quality score")

    # Test with noised
    noise_sim = ffd_sim.copy()
    noise_sim[:, 1] += 0.001 * rng.standard_normal(len(q))  # 5% noise

    quality_noisy_sim = FFQualityEvaluator.calc_ff_quality(noise_sim, ffd_sim)
    check.almost_equal(0.0, quality_noisy_sim, rel=1e-2, msg="Noisy simulation data should yield low quality score")

    # Experimental data: minimum shifted to q=0.3
    ffd_exp = np.column_stack([q, np.abs(np.sin(2 * np.pi * (q + 0.2))) * np.exp(-q)]).astype(float)

    quality = FFQualityEvaluator.calc_ff_quality(ffd_sim, ffd_exp)
    check.almost_equal(20, quality, rel=1e-2, msg="Shifted experimental data should yield +20 quality score")

    noise_exp = ffd_exp.copy()
    noise_exp[:, 1] += 0.001 * rng.standard_normal(len(q))  # 5% noise

    quality_noisy = FFQualityEvaluator.calc_ff_quality(noise_sim, noise_exp)
    check.almost_equal(
        20,
        quality_noisy,
        rel=1e-2,
        msg="Noisy simulation and experimental data should yield similar quality score to clean data",
    )


class _TestOPQualityBase:
    @pytest.fixture
    def mock_content(self):
        class MockLIPID:
            def __init__(slf, id):
                slf.id = id

            @property
            def mapping_dict(slf):
                return self.mapping_content[slf.id]

        contents = {k: MockLIPID(k) for k in self.mapping_content}
        yield contents

    @pytest.fixture
    def mock_explist(self):
        class MockExpList:
            def loc(slf, id):
                class Data:
                    data = self.exp_data[id]

                return Data()

        yield MockExpList()

    @pytest.fixture
    def mock_sim(self, mock_content):
        class MockSim(dict):
            def __init__(slf):
                super().__init__()
                slf["path"] = "mock/path"
                slf["EXPERIMENT"] = {"ORDERPARAMETER": {"DPPC": list(self.exp_data.keys())}}
                slf.op_data = self.sim_data
                slf.lipids = mock_content

            def membrane_composition(self, basis):
                return {"DPPC": 1}

        yield MockSim()


class TestOPQualityEvaluator1L1E(_TestOPQualityBase):
    @pytest.fixture(autouse=True)
    def setup_common(self):
        self.mapping_content = {
            "DPPC": {
                "M_C1_M": {
                    "ATOMNAME": "C1",
                    "FRAGMENT": "headgroup",
                },
                "M_C1H1_M": {
                    "ATOMNAME": "H1",
                    "FRAGMENT": "headgroup",
                },
                "M_C2_M": {
                    "ATOMNAME": "C2",
                    "FRAGMENT": "sn-1",
                },
                "M_C2H1_M": {
                    "ATOMNAME": "H2",
                    "FRAGMENT": "sn-1",
                },
            }
        }
        self.exp_data = {
            "exp1": {
                "DPPC": {
                    "M_C1_M M_C1H1_M": [-0.22, 0.02],
                    "M_C2_M M_C2H1_M": [-0.18, 0.01],
                },
            }
        }
        self.sim_data = {
            "DPPC": {
                "M_C1_M M_C1H1_M": [-0.21, 0.02, 0.02],
                "M_C2_M M_C2H1_M": [-0.17, 0.01, 0.01],
            }
        }

    def test_opqe_evaluate(self, mock_sim, mock_explist):
        from fairmd.lipids.quality import OPQualityEvaluator

        opq = OPQualityEvaluator(mock_sim, mock_explist)
        opq.evaluate_one()

        # test atomic qualities
        laq = opq.lipid_atomic_qualities
        check.is_instance(laq, dict)
        check.is_instance(laq["DPPC"], dict)
        check.equal(len(laq["DPPC"]), 1)
        check.equal(len(laq["DPPC"]["exp1"]), 2)
        check.almost_equal(laq["DPPC"]["exp1"].values(), [0.46, 0.65], abs=1e-2)
        print(laq)


class TestOPQualityEvaluator1L2E(_TestOPQualityBase):
    @pytest.fixture(autouse=True)
    def setup_common(self):
        self.mapping_content = {
            "DPPC": {
                "M_C1_M": {
                    "ATOMNAME": "C1",
                    "FRAGMENT": "headgroup",
                },
                "M_C1H1_M": {
                    "ATOMNAME": "H1",
                    "FRAGMENT": "headgroup",
                },
                "M_C2_M": {
                    "ATOMNAME": "C2",
                    "FRAGMENT": "sn-1",
                },
                "M_C2H1_M": {
                    "ATOMNAME": "H2",
                    "FRAGMENT": "sn-1",
                },
            },
        }
        self.exp_data = {
            "exp1": {
                "DPPC": {
                    "M_C1_M M_C1H1_M": [-0.22, 0.02],
                    "M_C2_M M_C2H1_M": [-0.18, 0.01],
                },
            },
            "exp2": {
                "DPPC": {
                    "M_C1_M M_C1H1_M": [-0.1, 0.02],
                    "M_C2_M M_C2H1_M": [-0.05, 0.01],
                },
            },
        }
        self.sim_data = {
            "DPPC": {
                "M_C1_M M_C1H1_M": [-0.21, 0.02, 0.02],
                "M_C2_M M_C2H1_M": [-0.17, 0.01, 0.01],
            },
        }

    def test_opqe_evaluate(self, mock_sim, mock_explist):
        from fairmd.lipids.quality import OPQualityEvaluator

        opq = OPQualityEvaluator(mock_sim, mock_explist)
        opq.evaluate_one()

        # test atomic qualities
        laq = opq.lipid_atomic_qualities
        check.is_instance(laq, dict)
        check.is_instance(laq["DPPC"], dict)
        check.equal(len(laq["DPPC"]), 2)
        check.equal(len(laq["DPPC"]["exp2"]), 2)
        print(laq)
