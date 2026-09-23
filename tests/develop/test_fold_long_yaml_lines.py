import importlib.util
import sys
from pathlib import Path

import pytest
import yaml

# developer/fold_long_yaml_lines.py is not part of the distributed package; see
# test_autocomplete_expsim_metadata.py for why these tests are marked `develop`.
pytestmark = [pytest.mark.develop, pytest.mark.nodata]

DEVELOPER_DIR = Path(__file__).resolve().parents[2] / "developer"

PROSE = (
    "Lipid dissolved in chloroform was dried under a stream of nitrogen and high vacuum overnight, "
    "and multilamellar liposomes were prepared by adding deuterium-depleted water."
)


def load_fold_module():
    """The script, loaded by path the way it runs: developer/ is on sys.path."""
    if str(DEVELOPER_DIR) not in sys.path:
        sys.path.insert(0, str(DEVELOPER_DIR))
    spec = importlib.util.spec_from_file_location("fold_long_yaml_lines", DEVELOPER_DIR / "fold_long_yaml_lines.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_a_long_value_is_folded_and_its_comment_kept():
    fold = load_fold_module()
    text = (
        "# hand-written header\n"
        f"SAMPLE_PROTOCOL: {PROSE}  # from the SI\n"
        "NMR:\n"
        f"  DETAILS: '{PROSE}'\n"
        "TEMPERATURE: 298\n"
    )
    out, left = fold.fold_file(text)
    lines = out.splitlines()

    assert yaml.safe_load(out) == yaml.safe_load(text)
    assert left == []
    assert all(len(line) <= fold.FOLD_WIDTH for line in lines)
    assert "SAMPLE_PROTOCOL: >-  # from the SI" in lines
    assert "  DETAILS: >-" in lines
    assert lines[0] == "# hand-written header"
    assert lines[-1] == "TEMPERATURE: 298"


def test_an_existing_folded_block_is_rewrapped():
    fold = load_fold_module()
    text = f"SAMPLE_PROTOCOL: >\n  {PROSE} {PROSE}\n\n  {PROSE}\n\nTEMPERATURE: 298\n"
    out, left = fold.fold_file(text)

    assert yaml.safe_load(out) == yaml.safe_load(text)
    assert left == []
    assert all(len(line) <= fold.FOLD_WIDTH for line in out.splitlines())
    assert out.startswith("SAMPLE_PROTOCOL: >\n")
    # The blank line before the next key is the author's, and stays.
    assert out.endswith("\n\nTEMPERATURE: 298\n")


def test_what_cannot_be_folded_is_left_and_reported():
    fold = load_fold_module()
    unbreakable = "(innerleaflet,76CHOL,30POPC,52POPE)(outerleaflet,76CHOL,52POPC,6POPE,2POPS,57PSM),13495SOL,52SOD"
    text = f"DETAILS: |\n  {PROSE}\nSYSTEM: {unbreakable}\nSHORT: fine\n"
    out, left = fold.fold_file(text)

    assert out == text
    assert left == [(1, "literal block"), (3, "no space to fold at")]
