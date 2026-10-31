"""
`test_misc` contains unit tests of auxiliary functions.

Test data is stored in `./ToyData/Simulations.1`

-------------------------------------------------------------------------------
NOTE: globally import of fairmd-lipids is **STRICTLY FORBIDDEN** because it
      breaks the substitution of global path folders
"""

import json

import numpy as np
import numpy.testing as npt
import pytest
import pytest_check as check

# run only on sim2 mocking data
pytestmark = [pytest.mark.sim1, pytest.mark.min]


def test_uname2element():
    """Test uname2element function."""
    from fairmd.lipids.auxiliary.mollib import uname2element

    check.equal(uname2element("M_C1_M"), "C")
    check.equal(uname2element("M_G1_M"), "C")
    check.equal(uname2element("M_C1N3_M"), "N")
    check.equal(uname2element("M_X_M"), "Dummy")
    check.equal(uname2element("M_D_M"), "Dummy")
    check.equal(uname2element("M_G1D_M"), "Dummy")
    check.equal(uname2element("M_G1O1D2_M"), "Dummy")

    with pytest.raises(KeyError):
        uname2element("UnknownElement")


@pytest.mark.parametrize(
    ("mappings", "composition", "expected"),
    [
        (
            {"DPPC": {"M_G1C4_M": {"ATOMNAME": "C4"}, "M_G1C5_M": {"ATOMNAME": "C5"}}},
            ("DPPC",),
            ("C5", ""),
        ),
        (
            {"DPPC": {"M_G3_M": {"ATOMNAME": "GLY3"}}},
            ("DPPC",),
            ("", "GLY3"),
        ),
        ({"DPPC": {"M_C1_M": {"ATOMNAME": "C1"}}}, ("DPPC",), ("", "")),
        (
            {
                "OTHER": {"M_G3_M": {"ATOMNAME": "OTHER_G3"}},
                "DPPC": {
                    "M_G3_M": {"ATOMNAME": "GLY3"},
                    "M_G1C4_M": {"ATOMNAME": "C4"},
                    "M_G1C5_M": {"ATOMNAME": "C5"},
                },
            },
            ("OTHER", "DPPC"),
            ("C5", "GLY3"),
        ),
    ],
)
def test_find_terminal_tail_and_g3_atoms(logger, mappings, composition, expected):
    """Test finding tail and glycerol carbon names for mapping scenarios."""
    from fairmd.lipids.auxiliary.mollib import find_terminal_tail_and_g3_atoms

    class MockLipid:
        def __init__(self, mapping_dict):
            self.mapping_dict = mapping_dict

    class MockSystem(dict):
        ID = 9999

        @property
        def content(self):
            return self["CONTENT"]

    content = {molecule: MockLipid(mapping) for molecule, mapping in mappings.items()}
    system = MockSystem(
        ID=9999,
        COMPOSITION={molecule: {"COUNT": 1} for molecule in composition},
        CONTENT=content,
    )

    assert find_terminal_tail_and_g3_atoms(system, logger) == expected


def test_get_tails():
    """Test Lipid.get_tails function."""
    from fairmd.lipids.auxiliary.mollib import get_tails_of_lipid
    from fairmd.lipids.molecules import Lipid

    # It counts tails for a phospholipid
    lipid = Lipid(name="DPPC")
    lipid.register_mapping("mappingDPPCberger.yaml")

    tails = get_tails_of_lipid(lipid)
    check.is_instance(tails, set)
    check.equal(len(tails), 2)
    check.is_in("sn-1", tails)
    check.is_in("sn-2", tails)

    # It counts tails for even cardiolipin
    lipid = Lipid(name="TOCL")
    lipid.register_mapping("mappingTOCLcharmm.yaml")

    tails = get_tails_of_lipid(lipid)
    check.equal(len(tails), 4)
    check.is_in("sn-1 1", tails)
    check.is_in("sn-2 1", tails)

    # It properly count tails for sphingolipids
    class MockLipid:
        mapping_dict: dict
        # obj

    mkl = MockLipid()
    mkl.mapping_dict = {
        "M_N1C1_M": {"FRAGMENT": "fa"},
        "M_C1_M": {"FRAGMENT": "sphingosine"},
        "M_P_M": {"FRAGMENT": "headgroup"},
    }
    tails = get_tails_of_lipid(mkl)
    check.equal(len(tails), 2)


def test_json_encoder(tmpdir):
    """Fixture for a mock experiment path with no data files."""
    from fairmd.lipids.auxiliary import CompactJSONEncoder

    exp_dir = tmpdir.mkdir("jsonenc")

    data = {
        "compact_object": {"first": "element", "second": 2},
        "compact_list": ["first", "second"],
        "long_list": [
            "this",
            "is",
            "a",
            "rather",
            "long\nlist",
            "and should be broken up because of its width",
        ],
        "non_ascii": "汉语",
        1: 2,
    }
    with open(exp_dir.join("test.yaml"), "w", encoding="utf-8") as fd:
        json.dump(data, fd, cls=CompactJSONEncoder, ensure_ascii=False)


def test_average_block():
    from fairmd.lipids.auxiliary import block_average_time_series

    # times 0..9, values = times
    t = np.arange(10, dtype=float)
    x = t.copy()
    arr = np.column_stack((t, x))

    out = block_average_time_series(arr, blocksize=2.0)
    # expected bins: [0,2), [2,4), [4,6), [6,8), [8,10]
    expected_times = np.array([1, 3, 5, 7, 9], dtype=float)
    expected_vals = np.array(
        [
            (0 + 1) / 2,
            (2 + 3) / 2,
            (4 + 5) / 2,
            (6 + 7) / 2,
            (8 + 9) / 2,
        ],
        dtype=float,
    )

    npt.assert_allclose(out[:, 0], expected_times)
    npt.assert_allclose(out[:, 1], expected_vals)

    arr = arr[:-1, :]
    out = block_average_time_series(arr, blocksize=2.0)


def test_progress(capsys):
    from fairmd.lipids import progress

    for i in progress(range(10)):
        pass
    captured = capsys.readouterr()
    check.is_in("100%", captured.out)
    check.is_not_in("100%", captured.err)
