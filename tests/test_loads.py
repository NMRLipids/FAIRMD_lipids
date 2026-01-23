"""
`test_loads` tests ONLY functions related to downloading files and/or resolving links.

NOTE: globally import of fairmd-lipids is **STRICTLY FORBIDDEN** because it
      breaks the substitution of global path folders
"""

import hashlib
import os
import requests
import responses
from responses import matchers
from unittest import mock
import sys
import time

import pytest
import pytest_check as check

# run only without mocking data
pytestmark = [pytest.mark.nodata, pytest.mark.min]


class TestDownloadWithProgressWithRetry:
    url = "https://example.org/file.bin"
    fname = "file.bin"

    @responses.activate
    def test_download_success(self, tmp_path):
        import fairmd.lipids.databankio as dio

        body = b"abcdef"
        dest = os.path.join(str(tmp_path), self.fname)

        responses.add(
            responses.GET,
            self.url,
            status=200,
            body=body,
            headers={"Content-Length": str(len(body))},
        )

        with check.raises(OSError) as e:
            dio.download_with_progress_with_retry("https://nonexistent.domain/file.bin", str(tmp_path))

        dio.download_with_progress_with_retry(self.url, dest)

        check.is_true(os.path.isfile(dest), "Success download must create a file")
        check.equal(open(dest, "br").read(), body, "Success download must get predefined content")

    @pytest.mark.parametrize(
        "maxsize, fsize, dsize",
        [
            (1000, 1500, 1000),
            (1000, 700, 700),
            (8192, 8192 * 3, 8192),  # test exactly chunk size
            (1000, 1000, 1000),
        ],
    )
    @responses.activate
    def test_download_dry_run(self, tmp_path, fsize, dsize, maxsize):
        import fairmd.lipids.databankio as dio

        dest = os.path.join(str(tmp_path), self.fname)
        body = b"x" * fsize

        responses.add(
            responses.GET,
            self.url,
            status=200,
            body=body,
            headers={"Content-Length": str(len(body))},
        )

        status = dio.download_with_progress_with_retry(
            self.url,
            dest,
            stop_after=maxsize,
        )

        check.is_true(os.path.isfile(dest), "Dry-run mode must create a file")
        check.equal(os.path.getsize(dest), dsize, "Stop-after mode must download not more than some number of bytes")

    @responses.activate
    def test_download_after_break(self, tmp_path):
        import fairmd.lipids.databankio as dio

        body_part1 = b"a" * dio._fmdl_chunk_size
        body_part2 = b"b" * dio._fmdl_chunk_size

        dest = os.path.join(str(tmp_path), self.fname)

        responses.add(
            responses.GET,
            self.url,
            status=200,
            body=body_part1 + body_part2,
            headers={"Content-Length": str(dio._fmdl_chunk_size * 2)},
        )
        responses.add(
            responses.GET,
            self.url,
            status=206,
            body=body_part2,
            headers={"Content-Length": str(dio._fmdl_chunk_size)},
            match=[matchers.header_matcher({"Range": "bytes=8192-"})],
        )

        def _ic_brok(chunk_size):
            yield body_part1
            raise requests.exceptions.ReadTimeout("timeout")

        with mock.patch(
            "requests.models.Response.iter_content",
            side_effect=_ic_brok,
        ):
            # first attempt, broken in the middle
            with check.raises(requests.exceptions.ReadTimeout):
                s = dio.download_with_progress_with_retry(self.url, dest)
            check.is_true(os.path.isfile(dest), "Partial download must create a file")
            check.equal(
                os.path.getsize(dest),
                dio._fmdl_chunk_size,
                "Partial download must have partial size",
            )

        # now resume
        dio.download_with_progress_with_retry(self.url, dest)

        check.is_true(os.path.isfile(dest), "Resumed download must create a file")
        check.equal(
            os.path.getsize(dest),
            dio._fmdl_chunk_size * 2,
            "Resumed download must have full size",
        )
        with open(dest, "br") as f:
            check.equal(
                f.read(),
                body_part1 + body_part2,
                "Resumed download must have full content",
            )


class TestDownloadResourceFromUri:
    url = "https://example.org/file.bin"
    fname = "file.bin"

    @responses.activate
    def test_download_success(self, tmp_path):
        import fairmd.lipids.databankio as dio

        body = b"abcdef"
        dest = os.path.join(str(tmp_path), self.fname)

        responses.add(
            responses.GET,
            self.url,
            status=200,
            body=body,
            headers={"Content-Length": str(len(body))},
        )

        status = dio.download_resource_from_uri(self.url, dest)

        check.equal(status, 0, "Success download must return zero")
        check.is_true(os.path.isfile(dest), "Success download must create a file")
        check.equal(open(dest, "br").read(), body, "Success download must get predefined content")

    @pytest.mark.parametrize(
        "fsize, dsize",
        [
            (8192 * 4, 0),
            (-8192 * 4, -8192 * 4),
        ],
    )
    @responses.activate
    def test_download_dry_run(self, tmp_path, fsize, dsize):
        import fairmd.lipids.databankio as dio

        dest = os.path.join(str(tmp_path), self.fname)
        body = b"x" * (dio.MAX_BYTES_DEFAULT + fsize)

        responses.add(
            responses.GET,
            self.url,
            status=200,
            body=body,
            headers={"Content-Length": str(len(body))},
        )

        status = dio.download_resource_from_uri(
            self.url,
            dest,
            max_bytes=True,
        )

        check.equal(status, 0, "Dry-run mode must work")
        check.is_true(os.path.isfile(dest), "Dry-run mode must create a file")
        check.equal(
            os.stat(dest).st_size,
            dio.MAX_BYTES_DEFAULT + dsize,
            "Dry-run mode must download not more than some number of bytes",
        )

    @responses.activate
    def test_resume_break_in_the_middle(self, tmp_path):
        import fairmd.lipids.databankio as dio

        body_part1 = b"a" * 8192
        body_part2 = b"b" * 8192

        dest = os.path.join(str(tmp_path), self.fname)
        responses.add(
            responses.GET,
            self.url,
            status=200,
            body=body_part1 + body_part2,
            headers={"Content-Length": str(16384)},
        )
        responses.add(
            responses.GET,
            self.url,
            status=206,
            body=body_part2,
            headers={"Content-Length": str(16384)},
            match=[matchers.header_matcher({"Range": "bytes=8192-"})],
        )

        def _ic_midbrok(*args, **kwargs):
            yield body_part1
            raise requests.exceptions.ReadTimeout("timeout")

        with mock.patch(
            "requests.models.Response.iter_content",
            side_effect=[_ic_midbrok(), iter([body_part2])],
        ):
            # now try again - should resume
            status = dio.download_resource_from_uri(self.url, dest, max_restarts=1)
            check.equal(status, 0, "Resume must work")
            check.is_true(os.path.isfile(dest), "Resume must create a file")
            check.is_false(os.path.isfile(dest + ".part"), "Resume must remove file.part")
            check.equal(
                os.path.getsize(dest),
                16384,
                "Resumed download must have full size",
            )

    @responses.activate
    def test_resume_nonmultiply(self, tmp_path):
        import fairmd.lipids.databankio as dio

        body_part1 = b"a" * 15000
        body_part2 = b"b" * 5000

        dest = os.path.join(str(tmp_path), self.fname)

        # first create a .part file with non-multiple of chunk size
        part_path = dest + ".part"
        with open(part_path, "wb") as f:
            f.write(body_part1)

        responses.add(  # for size-request mocking
            responses.GET,
            self.url,
            status=200,
            body=body_part1 + body_part2,
            headers={"Content-Length": str(20000)},
        )

        responses.add(
            responses.GET,
            self.url,
            status=206,
            body=body_part1[8192:] + body_part2,
            headers={"Content-Length": str(20000 - 8192)},
            match=[matchers.header_matcher({"Range": "bytes=8192-"})],
        )

        dio.download_resource_from_uri(self.url, dest)  # will truncate and continue
        with open(dest, "br") as fd:
            check.equal(fd.read(), body_part1 + body_part2, "Success download must get predefined content")

    @responses.activate
    def test_download_break_in_the_middle(self, tmp_path):
        import fairmd.lipids.databankio as dio

        body_part1 = b"a" * 8192
        body_part2 = b"b" * 8192

        def _ic_midbrok(*args, **kwargs):
            yield body_part1
            raise requests.exceptions.ReadTimeout("timeout")

        dest = os.path.join(str(tmp_path), self.fname)
        responses.add(
            responses.GET,
            self.url,
            status=200,
            body=body_part1 + body_part2,
            headers={"Content-Length": str(16384)},
        )

        with mock.patch(
            "requests.models.Response.iter_content",
            side_effect=[_ic_midbrok()],
        ):
            if os.path.exists(dest):
                os.remove(dest)
            with check.raises(ConnectionError):
                status = dio.download_resource_from_uri(self.url, dest)
            check.is_true(os.path.isfile(dest + ".part"), "Partial download must create a file.part")
            check.is_false(os.path.isfile(dest), "Partial download must not create a file.")


class TestGetFileSize:
    url = "https://example.org/file.bin"

    @responses.activate
    def test_no_content_length(self):
        import fairmd.lipids.databankio as dio

        responses.add(
            responses.GET,
            self.url,
            status=200,
            headers={},  # no Content-Length
        )

        size = dio._get_file_size_with_retry(self.url)

        assert size == 0

    @responses.activate
    def test_ok(self):
        import fairmd.lipids.databankio as dio

        responses.add(
            responses.GET,
            self.url,
            status=200,
            headers={"Content-Length": "1234"},
        )

        size = dio._get_file_size_with_retry(self.url)

        assert size == 1234
        assert len(responses.calls) == 1


class TestResolveFileUrl:
    def test_badDOI(self):
        import fairmd.lipids.databankio as dio

        time.sleep(5)
        # test if bad DOI fails
        print("Testing bad DOI resolution", file=sys.stderr)
        with check.raises(requests.exceptions.HTTPError) as e:
            dio.resolve_file_url("10.5281/zenodo.8435a", "a", validate_uri=True)
        check.is_in("404", str(e.value), "Bad zenodo ID should raise 404 error")
        # bad DOI doesn't fail if not to check
        print("Testing DOI resolution without validation", file=sys.stderr)
        check.equal(
            dio.resolve_file_url("10.5281/zenodo.8435a", "a.txt", validate_uri=False),
            "https://zenodo.org/records/8435a/files/a.txt",
        )
        time.sleep(5)
        # non-zenodo DOI fails
        with check.raises(NotImplementedError) as e:
            dio.resolve_file_url("10.1000/xyz123", "a.txt", validate_uri=False)
        check.is_in("Repository not validated", str(e.value))

    def test_goodDOI(self):
        import fairmd.lipids.databankio as dio

        time.sleep(5)
        # good DOI works properly
        assert (
            dio.resolve_file_url("10.5281/zenodo.8435138", "pope-md313rfz.tpr", validate_uri=True)
            == "https://zenodo.org/records/8435138/files/pope-md313rfz.tpr"
        )

    @pytest.mark.parametrize(
        "name, statuses, expected_exception",
        [
            ("transient 503 succeeds", [503, 503, 200], None),
            ("persistent 503 fails", [503] * 200, requests.exceptions.RetryError),
            ("403 fails immediately", [403], requests.exceptions.HTTPError),
        ],
    )
    @responses.activate
    def test_retry_logic(self, name, statuses, expected_exception):
        import fairmd.lipids.databankio as dio

        print(f"Testing resolve_doi_url with {name}", file=sys.stderr)
        url = "https://zenodo.org/api/records/8435138"

        for status in statuses:
            responses.add(responses.GET, url, status=status, body=b'{"files": [{"key": "a.txt"}]}')

        if expected_exception:
            with pytest.raises(expected_exception):
                dio.resolve_file_url("10.5281/zenodo.8435138", "a.txt", validate_uri=True)
        else:
            dio.resolve_file_url("10.5281/zenodo.8435138", "a.txt", validate_uri=True)

            assert len(responses.calls) == min(10, len(statuses))


class TestCalcFileSha1Hash:
    @staticmethod
    def _expected_sha1(data: bytes) -> str:
        return hashlib.sha1(data).hexdigest()

    @pytest.mark.parametrize(
        "data, step, expected_slice",
        [
            (b"hello world", 64, b"hello world"),
            (b"a" * 100, 10, b"a" * 10),
        ],
    )
    def test_one_block_behavior(self, tmp_path, data, step, expected_slice):
        from fairmd.lipids.databankio import calc_file_sha1_hash

        fi = os.path.join(str(tmp_path), "test.bin")
        with open(fi, "wb") as f:
            f.write(data)

        # one-block should be default
        result = calc_file_sha1_hash(fi, step=step)
        assert result == self._expected_sha1(expected_slice)

        result = calc_file_sha1_hash(fi, step=step, one_block=True)
        assert result == self._expected_sha1(expected_slice)

    @pytest.mark.parametrize(
        "data, step",
        [
            (b"abc" * 1000, 64),
            (b"hello", 1024),
        ],
    )
    def test_multi_block_reads_entire_file(self, tmp_path, data, step):
        from fairmd.lipids.databankio import calc_file_sha1_hash

        fi = os.path.join(str(tmp_path), "test.bin")
        with open(fi, "wb") as f:
            f.write(data)

        result = calc_file_sha1_hash(fi, step=step, one_block=False)

        assert result == self._expected_sha1(data)

    def test_empty_file(self, tmp_path):
        from fairmd.lipids.databankio import calc_file_sha1_hash

        fi = os.path.join(str(tmp_path), "test.bin")
        with open(fi, "wb") as f:
            pass

        with pytest.raises(ValueError):
            calc_file_sha1_hash(fi)


# TODO: create_simulation_directories
