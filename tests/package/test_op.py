"""
Test Order Parameter handling functionality.

Currently focused on testing opconvertor realted code because it is important for
OP data visualization.

Test data is stored in `./ToyData/Simulations.2`

-------------------------------------------------------------------------------
NOTE: globally import of fairmd-lipids is **STRICTLY FORBIDDEN** because it
      breaks the substitution of global path folders
"""

import logging
import MDAnalysis as mda
import numpy as np
import pytest
import pytest_check as check

# run only on sim2 mocking data
pytestmark = [pytest.mark.sim2, pytest.mark.min]


class TestBuildNiceOPdict:
    @pytest.fixture
    def systems(self):
        from fairmd.lipids.core import initialize_databank

        return initialize_databank()

    def test_work_empty(self):
        from fairmd.lipids.auxiliary.opconvertor import build_nice_OPdict
        from fairmd.lipids.molecules import Lipid

        lipid = Lipid("POPE")
        lipid.register_mapping()
        rdict = build_nice_OPdict({}, lipid)
        assert rdict == {}

    def test_build_fragmented(self, systems):
        from fairmd.lipids.auxiliary.opconvertor import build_nice_OPdict
        from fairmd.lipids.api import get_OP

        sys = systems.loc(281)
        opdata = get_OP(sys)
        rdict = build_nice_OPdict(opdata["POPC"], sys.lipids["POPC"])
        assert isinstance(rdict, dict)  # dict expected
        assert "sn-1" in rdict
        assert "sn-2" in rdict  # fragments at the top level

    def test_mockdicts(self):
        from fairmd.lipids.auxiliary.opconvertor import build_nice_OPdict
        from fairmd.lipids.molecules import Lipid

        lipid = Lipid("POPE")
        lipid.register_mapping()
        mock_opdata = {
            "M_G1C3_M M_G1C3H2_M": [0.1, 0.01, 0.001],
            "M_G1C3_M M_G1C3H1_M": [0.1, 0.01, 0.001],
            "M_G1C4_M M_G1C4H1_M": [0.2, 0.02, 0.002],
            "M_G1C4_M M_G1C4H2_M": [0.2, 0.02, 0.002],
            "M_G1C5_M M_G1C5H1_M": [0.3, 0.03, 0.003],
        }
        rdict = build_nice_OPdict(mock_opdata, lipid)
        check.is_in("sn-1", rdict)
        check.is_not_in("sn-2", rdict)
        check.equal(len(rdict["sn-1"]), 5)
        h_order = [int(x["H"]) for x in rdict["sn-1"]]
        check.is_true(h_order[:6] == [1, 2, 1, 2, 1], "sn-1 H ordering is not sorted")

    @pytest.mark.parametrize(
        ("values", "error_index", "expected"),
        [
            ([-0.1368, 0.00775834, 0.00048208], 2, 0.00048208),
            ([-0.1368, 0.00775834, 0.00048208], None, 0.00775834),
            ([-0.13, 0.02], None, 0.02),
            ([-0.13], None, None),
            ([-0.1368, 0.00775834], 2, None),
            ([-0.1368, 0.00775834, None], 2, None),
            ([-0.1368, 0.00775834, float("nan")], 2, None),
            ([-0.1368, 0.00775834, 0.0], 2, 0.0),
        ],
    )
    def test_error_selection(self, values, error_index, expected):
        from fairmd.lipids.auxiliary.opconvertor import build_nice_OPdict
        from fairmd.lipids.molecules import Lipid

        lipid = Lipid("POPC")
        lipid.register_mapping()
        kwargs = {} if error_index is None else {"error_index": error_index}
        result = build_nice_OPdict({"M_G1_M M_G1H1_M": values}, lipid, **kwargs)
        assert result["glycerol backbone"] == [
            {"C": "g1", "H": "1", "OP": values[0], "STD": expected},
        ]

    def test_cnames_pl(self, systems):
        from fairmd.lipids.auxiliary.opconvertor import build_nice_OPdict
        from fairmd.lipids.api import get_OP

        sys = systems.loc(281)
        opdata = get_OP(sys)
        rdict = build_nice_OPdict(opdata["POPC"], sys.lipids["POPC"])

        # C check numbers
        def has_c(cname: str, flist: dict) -> bool:
            return any(_c["C"] == cname for _c in flist)

        check.is_true(has_c("2", rdict["sn-1"]))
        check.is_true(has_c("16", rdict["sn-1"]))
        check.is_false(has_c("1", rdict["sn-1"]))
        check.is_false(has_c("17", rdict["sn-1"]))
        check.is_true(has_c("2", rdict["sn-2"]))
        check.is_true(has_c("18", rdict["sn-2"]))
        check.is_false(has_c("1", rdict["sn-2"]))
        check.is_false(has_c("19", rdict["sn-2"]))
        # H check names
        check.is_true(all(_c["H"] in ["1", "2", "3"] for _c in rdict["sn-1"]))
        check.is_true(all(_c["H"] in ["1", "2", "3"] for _c in rdict["sn-2"]))
        # check backbone
        check.is_true(has_c("g1", rdict["glycerol backbone"]))
        check.is_true(has_c("g2", rdict["glycerol backbone"]))
        check.is_true(has_c("g3", rdict["glycerol backbone"]))
        # check headgroup
        check.is_true(has_c("α", rdict["headgroup"]))
        # check ordering
        c_order = [int(x["C"]) for x in rdict["sn-1"]]
        check.is_true(c_order == sorted(c_order), "sn-1 C ordering is not sorted")
        h_order = [int(x["H"]) for x in rdict["sn-1"]]
        check.is_true(h_order[:6] == [1, 2, 1, 2, 1, 2], "sn-1 H ordering is not sorted")


def test_find_op_rotating_ethylene(rotating_ethylene) -> None:
    """Calculate C-H OPs for ethylene rotating in the YZ plane."""
    from fairmd.lipids.analib.databankop import find_OP

    mapping, universe = rotating_ethylene
    result = find_OP(mapping, universe, "ETH")
    check.equal(len(result), 4)
    check.is_true(all(np.isfinite(op.avg_std_stem[0]) for op in result))
    check.is_true(all(np.isfinite(op.avg_std_stem[1]) for op in result))
    check.is_true(all(np.isfinite(op.avg_std_stem[2]) for op in result))
    check.almost_equal(result[0].avg_std_stem[0], 0.53535093, rel=1e-4)
    check.almost_equal(result[1].avg_std_stem[0], 0.65017423, rel=1e-4)
    logger = logging.getLogger(__name__)
    logger.debug("C-H OPs for ethylene rotating in the YZ plane:")
    logger.debug("".join(f"\n{op.name}: {op.avg_std_stem[0]:.8f}" for op in result))


@pytest.fixture
def rotating_ethylene() -> tuple[dict, mda.Universe]:
    """Create a mock Universe of ethylene rotating in the YZ plane. Returns mapping + Universe"""
    ethylene = {
        "C1": (-0.67, 0.0),
        "H1": (-1.23, 0.92),
        "H2": (-1.23, -0.92),
        "C2": (0.67, 0.0),
        "H3": (1.23, 0.92),
        "H4": (1.23, -0.92),
    }
    mapping = {
        "M_C1_M": {"ATOMNAME": "C1"},
        "M_C1H1_M": {"ATOMNAME": "H1"},
        "M_C1H2_M": {"ATOMNAME": "H2"},
        "M_C2_M": {"ATOMNAME": "C2"},
        "M_C2H1_M": {"ATOMNAME": "H3"},
        "M_C2H2_M": {"ATOMNAME": "H4"},
    }
    names = list(ethylene)
    universe = mda.Universe.empty(
        len(names),
        n_residues=1,
        atom_resindex=np.zeros(len(names), dtype=int),
        trajectory=True,
    )
    universe.add_TopologyAttr("name", names)
    universe.add_TopologyAttr("resname", ["ETH"])
    universe.add_TopologyAttr("resid", [1])
    angles = np.deg2rad(0.05 * np.arange(100))
    positions = np.zeros((100, len(names), 3), dtype=np.float32)
    for frame, angle in enumerate(angles):
        rotation = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
        positions[frame, :, 1:3] = np.array([rotation @ ethylene[name] for name in names])
    universe.load_new(positions, order="fac")
    return mapping, universe


def test_notebook_error_bars_use_simulation_sem_and_experimental_uncertainty():
    matplotlib = pytest.importorskip("matplotlib")
    pytest.importorskip("IPython")

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from fairmd.lipids.ipylib.plotop import _plot_glyhead_united
    from fairmd.lipids.molecules import Lipid

    lipid = Lipid("POPC")
    lipid.register_mapping()
    simulation = {"M_G1_M M_G1H1_M": [-0.1368, 0.00775834, 0.00048208]}
    experiment = {"M_G1_M M_G1H1_M": [-0.13, 0.02]}
    figures = _plot_glyhead_united(simulation, experiment, lipid)
    try:
        axis = figures["head"].axes[0]
        sim_bar, exp_bar = axis.collections
        np.testing.assert_allclose(sim_bar.get_segments()[0][:, 1], [-0.13728208, -0.13631792])
        np.testing.assert_allclose(exp_bar.get_segments()[0][:, 1], [-0.15, -0.11])
    finally:
        for figure in figures.values():
            plt.close(figure)
