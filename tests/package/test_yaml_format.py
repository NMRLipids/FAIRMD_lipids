"""
`test_yaml_format` tests ONLY YAML public formatting function..

NOTE: globally import of fairmd-lipids is **STRICTLY FORBIDDEN** because it
      breaks the substitution of global path folders
"""

import yaml

import pytest

pytestmark = [pytest.mark.nodata, pytest.mark.min]

PROSE = (
    "Lipid dissolved in chloroform was dried under a stream of nitrogen and high vacuum overnight, "
    "and multilamellar liposomes were prepared by adding deuterium-depleted water."
)


def test_a_long_value_is_folded_and_its_value_preserved():
    from fairmd.lipids.auxiliary import encode_canonical_yaml

    data = {
        "SAMPLE_PROTOCOL": PROSE,
        "NMR": {"DETAILS": PROSE},
        "TEMPERATURE": 298,
    }

    out = encode_canonical_yaml(data)
    lines = out.splitlines()

    assert yaml.safe_load(out) == data
    assert all(len(line) <= 100 for line in lines)
    assert "SAMPLE_PROTOCOL: >-" in lines
    assert "  DETAILS: >-" in lines
    assert lines[-1] == "TEMPERATURE: 298"


def test_existing_folded_block_is_canonicalized():
    from fairmd.lipids.auxiliary import encode_canonical_yaml

    data = {"SAMPLE_PROTOCOL": f"{PROSE} {PROSE}\n\n{PROSE}", "TEMPERATURE": 298}

    out = encode_canonical_yaml(data)

    assert yaml.safe_load(out) == data
    assert out.endswith("TEMPERATURE: 298\n")


def test_unbreakable_value_is_serialized_without_changing_it():
    from fairmd.lipids.auxiliary import encode_canonical_yaml

    unbreakable = "(innerleaflet,76CHOL,30POPC,52POPE)(outerleaflet,76CHOL,52POPC,6POPE,2POPS,57PSM),13495SOL,52SOD"
    data = {"DETAILS": PROSE, "SYSTEM": unbreakable, "SHORT": "fine"}

    out = encode_canonical_yaml(data)

    assert yaml.safe_load(out) == data
