"""
`test_loads` tests ONLY functions related to downloading files and/or resolving links.

NOTE: globally import of fairmd-lipids is **STRICTLY FORBIDDEN** because it
      breaks the substitution of global path folders
"""

import os
import stat
import sys
from urllib.error import HTTPError, URLError

import pytest
import pytest_check as check

# run only without mocking data
pytestmark = [pytest.mark.nodata, pytest.mark.min]


class TestDownloadResourceFromUri:
    TESTFILENAME = "t.tpr"

    def test_justdl__download_resource_from_uri(self, monkeypatch, tmp_path):
        monkeypatch.chdir(tmp_path)
        import fairmd.lipids.databankio as dio

        if os.path.exists(self.TESTFILENAME):
            os.remove(self.TESTFILENAME)  # just for sure
        # first-time download
        assert (
            dio.download_resource_from_uri(
                "https://zenodo.org/records/8435138/files/pope-md313rfz.tpr",
                "./" + self.TESTFILENAME,
            )
            == 0
        )
        # repeat download
        check.equal(
            dio.download_resource_from_uri(
                "https://zenodo.org/records/8435138/files/pope-md313rfz.tpr",
                "./" + self.TESTFILENAME,
            ),
            1,
        )
        os.remove(self.TESTFILENAME)

    def test_corrupted__download_resource_from_uri(self, monkeypatch, tmp_path):
        monkeypatch.chdir(tmp_path)
        import fairmd.lipids.databankio as dio

        # redownload corrupted file
        with open(self.TESTFILENAME, "w") as f:
            f.write("BABABA")
        old_size = os.path.getsize(self.TESTFILENAME)
        assert (
            dio.download_resource_from_uri(
                "https://zenodo.org/records/8435138/files/pope-md313rfz.tpr",
                "./" + self.TESTFILENAME,
            )
            == 2
        )
        # check filesize
        check.greater(
            os.path.getsize(self.TESTFILENAME),
            old_size,
            msg="File was not redownloaded despite size mismatch!",
        )
        os.remove(self.TESTFILENAME)

    @pytest.mark.skipif(sys.platform == "darwin", reason="fails on GitHub Actions")
    def test_errs__download_resource_from_uri(self, monkeypatch, tmp_path):
        monkeypatch.chdir(tmp_path)
        d = tmp_path / "sub"
        d.mkdir()

        import fairmd.lipids.databankio as dio

        # put directory instead of filename
        with pytest.raises(IsADirectoryError) as _:
            dio.download_resource_from_uri("https://zenodo.org/records/8435138/files/pope-md313rfz.tpr", "sub")

        if os.name != "nt":
            # TODO: fix on Windows
            d.chmod(stat.S_IREAD | stat.S_IRGRP | stat.S_IROTH)
            fn = os.path.join("sub", "no-rights.tpr")

            # ask to write to file which you don't have an access
            with pytest.raises(PermissionError) as _:
                dio.download_resource_from_uri("https://zenodo.org/records/8435138/files/pope-md313rfz.tpr", fn)


# resolve_doi_url
class TestResolveDoiUrl:
    def test_badDOI__resolve_doi_url(self):
        import fairmd.lipids.databankio as dio

        # test if bad DOI fails
        with pytest.raises(HTTPError, match="404") as _:
            dio.resolve_doi_url("10.5281/zenodo.8435a", True)
        # bad DOI doesn't fail if not to check
        assert dio.resolve_doi_url("10.5281/zenodo.8435a", False) == "https://doi.org/10.5281/zenodo.8435a"

    def test_goodDOI__resolve_doi_url(self):
        import fairmd.lipids.databankio as dio

        # good DOI works properly
        assert dio.resolve_doi_url("10.5281/zenodo.8435138", True) == "https://doi.org/10.5281/zenodo.8435138"

    @staticmethod
    def _create_mock_success_response(mocker):
        mock_response = mocker.MagicMock()
        mock_response.read.return_value = b"mocked content"
        mock_response.getheader.return_value = "100"  # Content-Length
        mock_response.code = 200
        mock_response.reason = "OK"
        return mock_response

    @pytest.mark.parametrize(
        "name, side_effects_func, expected_exception, expected_call_count",
        [
            (
                "transient URLError succeeds",
                lambda m: [URLError("err1"), URLError("err2"), TestResolveDoiUrl._create_mock_success_response(m)],
                None,
                3,
            ),
            (
                "persistent URLError fails",
                lambda m: URLError("persistent error"),
                ConnectionError,
                5,
            ),
            (
                "non-retriable 403 HTTPError fails immediately",
                lambda m: HTTPError("url", 403, "Forbidden", {}, None),
                HTTPError,
                1,
            ),
            (
                "retriable 503 HTTPError succeeds",
                lambda m: [
                    HTTPError("url", 503, "Service Unavailable", {}, None),
                    HTTPError("url", 503, "Service Unavailable", {}, None),
                    TestResolveDoiUrl._create_mock_success_response(m),
                ],
                None,
                3,
            ),
            (
                "persistent 503 HTTPError fails",
                lambda m: HTTPError("url", 503, "Service Unavailable", {}, None),
                ConnectionError,
                5,
            ),
        ],
    )
    def test_retry_logic__resolve_doi_url(
        self,
        name,
        side_effects_func,
        expected_exception,
        expected_call_count,
        mocker,
    ):
        import fairmd.lipids.databankio as dio

        mocker.patch("time.sleep", return_value=None)

        side_effects = side_effects_func(mocker)

        mock_urlopen = mocker.patch(
            "fairmd.lipids.databankio.urllib.request.urlopen",
            side_effect=side_effects,
        )

        if expected_exception:
            with pytest.raises(expected_exception) as excinfo:
                dio.resolve_doi_url("10.5281/zenodo.8435138", True)
            if expected_exception == HTTPError:
                actual_http_error = side_effects[0] if isinstance(side_effects, list) else side_effects
                assert excinfo.value.code == actual_http_error.code
        else:
            dio.resolve_doi_url("10.5281/zenodo.8435138", True)

        assert mock_urlopen.call_count == expected_call_count, f"Test '{name}' failed on call count."

    # resolve_download_file_url

    def test__resolve_download_file_url(self):
        pass
