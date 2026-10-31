"""
Test module test_formfactor.py

Contains tests of functions designed for processing of form-factor curve.

Test data is stored in `./ToyData/Simulations.2`

-------------------------------------------------------------------------------
NOTE: globally import of fairmd-lipids is **STRICTLY FORBIDDEN** because it
      breaks the substitution of global path folders
"""

import numpy as np

import pytest
import pytest_check as check

pytestmark = [pytest.mark.nodata, pytest.mark.min]


def test_get_mins_from_ffdata():
    from fairmd.lipids.analib.formfactor import get_mins_from_ffdata

    synt_data = np.zeros((1000, 3))
    synt_data[:, 0] = np.linspace(0, 1, 1000)
    synt_data[:, 1] = np.abs(np.sin(synt_data[:, 0] * 8))
    # sin(x*8) = 0
    # x1 = pi/8
    # x2 = pi/4
    p = get_mins_from_ffdata(synt_data)
    check.almost_equal(p[0], np.pi / 8, abs=1e-3)
    check.almost_equal(p[1], np.pi / 4, abs=1e-3)


def test_get_mins_from_ffdata_rejects_data_with_too_few_points():
    """A Savitzky-Golay quadratic cannot be fitted to the resulting window."""
    from fairmd.lipids.analib.formfactor import get_mins_from_ffdata

    ffdata = np.array([[0.0, 1.0], [0.1, 0.5], [0.2, 1.0]])
    with pytest.raises(ValueError, match="Problems running Savitsky-Golay"):
        get_mins_from_ffdata(ffdata)


def test_calc_ff_scaling_distance():
    """
    Test when simulation values are twice the experimental values.

    The least-squares scaling coefficient is two and the residual is zero.
    """
    from fairmd.lipids.analib.formfactor import calc_ff_scaling_distance

    q = np.linspace(0.1, 0.4, 4)
    exp = np.column_stack((q, [1.0, 2.0, 3.0, 4.0], [0.5] * 4))
    sim = np.column_stack((q, [2.0, 4.0, 6.0, 8.0]))

    scaling, chi = calc_ff_scaling_distance(exp, sim)

    check.almost_equal(scaling, 2.0)
    check.almost_equal(chi, 0.0, abs=1e-12)


def test_estimate_error_of_min():
    from fairmd.lipids.analib.formfactor import calc_minpos_with_error, get_mins_from_ffdata

    # synthetic data
    synt_data = np.zeros((1000, 3))
    synt_data[:, 0] = np.linspace(0, 1, 1000)
    synt_data[:, 1] = np.abs(np.sin(synt_data[:, 0] * 8))
    synt_data[:, 2] = 0.1
    # Case 0: no error
    pos_noerr, err_noerr = calc_minpos_with_error(synt_data[:, :2])
    # this algorithm gives slightly different value of min. But it must be close
    check.almost_equal(
        pos_noerr,
        get_mins_from_ffdata(synt_data[:, :2])[0],
        abs=5e-4,
        msg="Minimum precise position is far from the legacy algo",
    )
    # no error is handled as constant error 0.1
    pos_err01, err_err01 = calc_minpos_with_error(synt_data)
    check.almost_equal(
        err_noerr,
        err_err01,
        abs=1e-3,
        msg="Error of minimum with no error is not estimated as err=0.1 [const]",
    )
    # Case 1: 2 constant errors
    pos_comp, err_comp = calc_minpos_with_error(synt_data)
    synt_data[:, 2] = 0.2
    pos_berr, err_berr = calc_minpos_with_error(synt_data)
    check.less(
        err_comp,
        err_berr,
        msg="Error of minimum with constant error 0.2 must be > than with err 0.1",
    )
    check.almost_equal(
        pos_comp,
        pos_berr,
        abs=1e-5,
        msg="Minimum position with constant error 0.2 must be close to minimum position with constant error 0.1",
    )
    # Case 2: noise
    synt_data[:, 1] += 0.05 * np.random.rand(1000) - 0.025
    pos_noised, err_noised = calc_minpos_with_error(synt_data)
    check.less(
        err_comp,
        err_noised,
        msg="Error of minimum with noise is not greater than error of minimum with constant error",
    )
    check.almost_equal(
        pos_comp,
        pos_noised,
        abs=1e-3,
        msg="Minimum position with and w/o noise should be similar",
    )
