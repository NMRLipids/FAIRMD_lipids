"""
Perform integration testing of adding-simulation functionality.

NOTE: globally import of fairmd-lipids is **STRICTLY FORBIDDEN** because it
      breaks the substitution of global path folders
"""

import os
import shutil
import subprocess
import time
from tempfile import TemporaryDirectory

import pytest
from pytest_check import check

# run only without mocking data
pytestmark = pytest.mark.adddata


@pytest.fixture(scope="module")
def tmp_work_dir():
    with TemporaryDirectory(prefix="dbtestWD_", dir=os.path.dirname(__file__)) as wdir:
        print(f"Will use following directory for loadings: {wdir}")
        yield wdir


class TestAddData:
    @classmethod
    def setup_class(cls):
        from fairmd import lipids

        if os.path.isfile(os.path.join(lipids.FMDL_DATA_PATH, ".notest")):
            pytest.exit("Test are corrupted. I see '.notest' file in the data folder.")
        cls.exe = os.path.join(
            os.path.dirname(os.path.dirname(__file__)), "src", "fairmd", "lipids", "bin", "add_simulation.py"
        )
        cls.out_dir = lipids.FMDL_SIMU_PATH
        os.mkdir(cls.out_dir)

    @classmethod
    def teardown_class(cls):
        if os.path.exists(cls.out_dir):
            shutil.rmtree(cls.out_dir, ignore_errors=True)

    def test_add_data_h(self):
        """
        Testing `add_simulation.py -h` behavior
        """
        result = subprocess.run(
            [
                self.exe,
                "-h",
            ],
            check=False,
            capture_output=True,
            text=True,
        )
        print(result.stdout)
        assert "--dry-run", "Expected --dry-run option in help output"
        assert result.returncode == 0

    @pytest.mark.parametrize("infofn, debug", [("info566.yaml", False), ("info566.yaml", True)])
    def test_add_data_addgood(self, infofn, debug, tmp_work_dir, capsys, request):
        """
        Testing `add_simulation.py -f <filename> -w <dirname> -o <dirname>` which should
        end correctly
        """
        fn = os.path.join(os.path.dirname(__file__), "ToyData", "info", infofn)
        run_list = [
            self.exe,
            "-f",
            fn,
            "-w",
            tmp_work_dir,
        ]
        if debug:
            run_list.append("-d")
        result = subprocess.run(run_list, check=False, capture_output=True, text=True)
        if request.config.getoption("capture") == "no":
            # do not disable capturing if pytest is run with `-s` option
            with capsys.disabled():
                print("stderr:")
                print(result.stderr)
                print("stdout:")
                print(result.stdout)
        time.sleep(1)
        if debug:
            "[DEBUG]" in result.stderr
            check.is_in("[DEBUG]", result.stderr, msg="Expected [DEBUG] in stderr when debug mode is on")
        else:
            check.is_not_in("[DEBUG]", result.stderr, msg="Expected no [DEBUG] in stderr when debug mode is off")
        check.is_in("[INFO]", result.stderr, msg="Expected [INFO] in stderr")
        self._check_new_readme(capsys)
        check.equal(result.returncode, 0)
        TestAddData.teardown_class()  # clean up to run the same test again

    @pytest.mark.parametrize("infofn", ["info566_uf.yaml"])
    def test_add_data_fail(self, infofn, tmp_work_dir, capsys):
        fn = os.path.join(os.path.dirname(__file__), "ToyData", "info", infofn)
        result = subprocess.run(
            [
                self.exe,
                "-f",
                fn,
                "-w",
                tmp_work_dir,
                "-d",
            ],
            check=False,
            capture_output=True,
            text=True,
        )
        print(result.stderr)
        assert "[ERROR]" in result.stderr, "Expected an error message in stderr"
        assert result.returncode == 1

    def _check_new_readme(self, capsys):
        from fairmd.lipids.core import initialize_databank

        ss = initialize_databank()
        captured = capsys.readouterr()
        if "!!README LOAD ERROR!!" in captured.err:
            pytest.fail("There was an error in the captured output: " + captured.err)
        for s in ss:
            print(s["DOI"])
        check.is_true(len(ss) > 0, msg="No systems found in the databank after adding data.")
