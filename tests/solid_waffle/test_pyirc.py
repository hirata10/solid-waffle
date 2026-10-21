"""Tests for features in pyirc."""


import numpy as np
import pytest
from solid_waffle.pyirc import (
    IndexDictionary,
    gain_alphabetacorr,
    gain_alphacorr,
    pyIRC_percentile,
    slidemed_percentile,
)


def test_indexerr():
    """Checks index raises an error if you add things in the wrong order."""

    d = IndexDictionary(0)
    d.addhnl(4)
    with pytest.raises(ValueError):
        d.addbfe(2)


def test_percentile_cont():
    """Checks the behavior of the percentile function if discrete option is off."""

    arr = np.linspace(0, 99, 100)
    mask = np.ones(100, dtype=bool)
    x = pyIRC_percentile(arr, mask, 75, disc=False)
    assert 74.24 < x < 74.26


def test_slidemed():
    """Test for sliding median method: get m where percentile(mx-y, target)=0."""

    theta = np.linspace(0.0, 2.0 * np.pi, 1001)[:-1]
    x = 3 + 2 * np.cos(theta)
    y = 0.01 * (2.5 * x + np.sin(theta))
    vals = []
    for p in [25, 50, 75]:
        vals.append(slidemed_percentile(x, y, p))
    assert np.allclose(vals, [0.022328464759931173, 0.025, 0.02767153524006883])

    theta = np.linspace(0.0, 2.0 * np.pi, 1001)[:-1]
    x = 3 + 3.5 * np.cos(theta)
    y = 0.01 * (2.5 * x + np.sin(theta))
    vals = []
    for p in [25, 50, 75]:
        vals.append(slidemed_percentile(x, y, p))
    assert np.allclose(vals, [0.02083919936463357, 0.025, 0.029160800635366434])


def test_return_gain_alphacorr():
    """Test failure of gain_alpha_corr."""

    assert len(gain_alphacorr(1.0, 100.0, 100.0, 100.0)) == 0
    assert len(gain_alphabetacorr(1.0e-6, 0.0, 0.0, 1.0e3, 0.0, [1, 2, 3, 4])) == 0
