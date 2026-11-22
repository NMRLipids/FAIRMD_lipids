import importlib.util
import os
import subprocess
import sys
from pathlib import Path

import pytest

# Exercises developer/reanchor_eof_hunks.py, which lives in the `developer/`
# folder and is not part of the distributed package.
pytestmark = [pytest.mark.develop, pytest.mark.nodata]

DEVELOPER_DIR = Path(__file__).resolve().parents[2] / "developer"


def load_module():
    module_path = DEVELOPER_DIR / "reanchor_eof_hunks.py"
    spec = importlib.util.spec_from_file_location("reanchor_eof_hunks", module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


mod = load_module()


def git(repo, *args):
    env = {
        "PATH": os.defpath,
        "HOME": str(repo / ".git-home"),
        "LANG": "C.UTF-8",
    }
    return subprocess.run(
        ["git", "-c", "user.email=t@t", "-c", "user.name=t", *args],
        cwd=repo,
        env=env,
        capture_output=True,
        text=True,
        check=True,
    ).stdout


@pytest.fixture
def repo(tmp_path):
    """A one-file git repository whose file is committed and can then be edited."""
    git(tmp_path, "init", "-q", ".")
    (tmp_path / "README.yaml").write_text("A: 1\nB: 2\nC: 3\n", encoding="utf-8")
    git(tmp_path, "add", "-A")
    git(tmp_path, "commit", "-qm", "base")
    return tmp_path


def anchors(diff):
    """The ``@@`` headers, plus every changed line, as reviewdog would read them."""
    return [
        line
        for line in diff.splitlines()
        if line[:1] in "@-+" and not line.startswith("+++") and not line.startswith("---")
    ]


# ---------------------------------------------------------------------------
# What the filter is for
# ---------------------------------------------------------------------------


def test_an_append_becomes_a_replacement_of_the_last_line(repo):
    """The case the autocomplete produces: a block added below the final line.

    Unfiltered this is ``@@ -3,0 +4,2 @@``, which reviewdog anchors one line past
    the end of a three-line file. Afterwards the last line is deleted and re-added
    with the new content behind it, so the anchor is line 3, which exists.
    """
    (repo / "README.yaml").write_text("A: 1\nB: 2\nC: 3\nD: 4\nE: 5\n", encoding="utf-8")
    out = mod.reanchor(git(repo, "diff"))
    assert anchors(out) == ["@@ -1,3 +1,5 @@", "-C: 3", "+C: 3", "+D: 4", "+E: 5"]


def test_the_rewritten_diff_still_applies_cleanly(repo):
    """The suggestion has to produce the file the script wrote, byte for byte."""
    written = "A: 1\nB: 2\nC: 3\nD: 4\n"
    (repo / "README.yaml").write_text(written, encoding="utf-8")
    patch = repo / "re.diff"
    patch.write_text(mod.reanchor(git(repo, "diff")), encoding="utf-8")
    git(repo, "checkout", "--", "README.yaml")
    git(repo, "apply", str(patch))
    assert (repo / "README.yaml").read_text(encoding="utf-8") == written


def test_hunk_line_counts_stay_correct(repo):
    """git apply is strict about ``@@ -a,b +c,d @@``, so the counts must still add up.

    Turning one context line into a ``-``/``+`` pair spends one old line and one
    new line exactly as the context line did, which is why the header is left
    alone.
    """
    (repo / "README.yaml").write_text("A: 1\nB: 2\nC: 3\nD: 4\n", encoding="utf-8")
    before = git(repo, "diff")
    after = mod.reanchor(before)
    assert [x for x in before.splitlines() if x.startswith("@@")] == [
        x for x in after.splitlines() if x.startswith("@@")
    ]


# ---------------------------------------------------------------------------
# What it must leave alone
# ---------------------------------------------------------------------------


def test_an_insertion_in_the_middle_is_left_alone(repo):
    """Only an append needs re-anchoring; anywhere else the anchor already exists."""
    (repo / "README.yaml").write_text("A: 1\nB: 2\nX: 9\nC: 3\n", encoding="utf-8")
    diff = git(repo, "diff")
    assert mod.reanchor(diff) == diff


def test_a_rewritten_last_line_is_left_alone(repo):
    """Trailing additions that follow a deletion are already on a real line."""
    (repo / "README.yaml").write_text("A: 1\nB: 2\nC: 9\nD: 4\n", encoding="utf-8")
    diff = git(repo, "diff")
    assert mod.reanchor(diff) == diff
    assert "-C: 3" in diff


def test_a_missing_final_newline_is_left_alone(repo):
    """A ``\\ No newline`` marker stops the scan rather than being rewritten.

    git already expresses this as a change to the last line, so reviewdog anchors
    it correctly, and splicing around the marker would only corrupt the diff.
    """
    (repo / "README.yaml").write_text("A: 1\nB: 2\nC: 3\nD: 4", encoding="utf-8")
    diff = git(repo, "diff")
    assert "\\ No newline at end of file" in diff
    assert mod.reanchor(diff) == diff


def test_only_the_last_hunk_of_a_file_is_considered(repo):
    """An earlier hunk ending in additions is followed by context, not by the end."""
    (repo / "README.yaml").write_text("A: 1\n" + "pad\n" * 10 + "B: 2\nC: 3\n", encoding="utf-8")
    out = mod.reanchor(git(repo, "diff"))
    # The padding hunk is untouched; only the file's final hunk may be rewritten.
    assert out.count("@@") == git(repo, "diff").count("@@")
    assert "-C: 3" not in out


def test_each_file_is_re_anchored_independently(repo):
    """A multi-file diff must not let one file's last hunk leak into the next."""
    (repo / "other.yaml").write_text("P: 1\nQ: 2\n", encoding="utf-8")
    git(repo, "add", "-A")
    git(repo, "commit", "-qm", "two files")
    (repo / "README.yaml").write_text("A: 1\nB: 2\nC: 3\nD: 4\n", encoding="utf-8")
    (repo / "other.yaml").write_text("P: 1\nQ: 2\nR: 3\n", encoding="utf-8")
    out = mod.reanchor(git(repo, "diff"))
    assert "-C: 3" in out and "+C: 3" in out
    assert "-Q: 2" in out and "+Q: 2" in out


def test_an_empty_diff_survives(repo):
    assert mod.reanchor("") == ""


# ---------------------------------------------------------------------------
# The command line
# ---------------------------------------------------------------------------


def test_it_works_as_a_stdin_to_stdout_filter(repo):
    """The workflow pipes git diff through it, so the CLI shape matters."""
    (repo / "README.yaml").write_text("A: 1\nB: 2\nC: 3\nD: 4\n", encoding="utf-8")
    result = subprocess.run(
        [sys.executable, str(DEVELOPER_DIR / "reanchor_eof_hunks.py")],
        input=git(repo, "diff"),
        capture_output=True,
        text=True,
        check=True,
    )
    assert "-C: 3" in result.stdout
