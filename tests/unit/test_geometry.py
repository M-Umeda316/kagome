"""Tests for geometry utilities (orthorhombic PBC)."""
import numpy as np
import pytest

from kagome.geometry import minimum_image, wrap_positions


class TestMinimumImage:

    def test_no_cell_passthrough(self):
        r = np.array([5.0, 0.0, 0.0])
        assert np.array_equal(minimum_image(r, None), r)

    def test_orthorhombic_wraps(self):
        cell = np.diag([10.0, 10.0, 10.0])
        r = np.array([7.0, -6.0, 3.0])
        result = minimum_image(r, cell)
        assert result == pytest.approx([-3.0, 4.0, 3.0])

    def test_triclinic_raises(self):
        cell = np.array([[10.0, 1.0, 0.0],
                         [0.0, 10.0, 0.0],
                         [0.0, 0.0, 10.0]])
        with pytest.raises(ValueError, match='orthorhombic'):
            minimum_image(np.array([1.0, 0.0, 0.0]), cell)


class TestWrapPositions:

    def test_no_cell_noop(self):
        pos = np.array([[5.0, 0.0, 0.0]])
        wrap_positions(pos, None)
        assert pos[0, 0] == 5.0

    def test_orthorhombic_wraps_inplace(self):
        cell = np.diag([10.0, 10.0, 10.0])
        pos = np.array([[12.0, -1.0, 5.0]])
        wrap_positions(pos, cell)
        assert pos[0] == pytest.approx([2.0, 9.0, 5.0])

    def test_triclinic_raises(self):
        cell = np.array([[10.0, 0.0, 0.0],
                         [0.0, 10.0, 2.0],
                         [0.0, 0.0, 10.0]])
        pos = np.array([[1.0, 2.0, 3.0]])
        with pytest.raises(ValueError, match='orthorhombic'):
            wrap_positions(pos, cell)
