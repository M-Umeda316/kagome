"""Tests for geometry utilities (orthorhombic PBC)."""
import numpy as np
import pytest

from kagome.geometry import (
    minimum_image, wrap_positions,
    validated_box, minimum_image_fast, wrap_positions_fast,
)


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


class TestValidatedBox:

    def test_returns_diagonal(self):
        cell = np.diag([8.0, 9.0, 10.0])
        box = validated_box(cell)
        assert box == pytest.approx([8.0, 9.0, 10.0])

    def test_none_passthrough(self):
        assert validated_box(None) is None

    def test_triclinic_raises(self):
        cell = np.array([[10.0, 1.0, 0.0],
                         [0.0, 10.0, 0.0],
                         [0.0, 0.0, 10.0]])
        with pytest.raises(ValueError, match='orthorhombic'):
            validated_box(cell)


class TestMinimumImageFast:

    def test_matches_original(self):
        cell = np.diag([10.0, 10.0, 10.0])
        box = validated_box(cell)
        r = np.array([7.0, -6.0, 3.0])
        expected = minimum_image(r, cell)
        result = minimum_image_fast(r, box)
        np.testing.assert_allclose(result, expected)

    def test_none_box_passthrough(self):
        r = np.array([5.0, 0.0, 0.0])
        assert np.array_equal(minimum_image_fast(r, None), r)


class TestWrapPositionsFast:

    def test_matches_original(self):
        cell = np.diag([10.0, 10.0, 10.0])
        box = validated_box(cell)
        pos_orig = np.array([[12.0, -1.0, 5.0]])
        pos_fast = pos_orig.copy()
        wrap_positions(pos_orig, cell)
        wrap_positions_fast(pos_fast, box)
        np.testing.assert_allclose(pos_fast, pos_orig)

    def test_none_box_noop(self):
        pos = np.array([[5.0, 0.0, 0.0]])
        wrap_positions_fast(pos, None)
        assert pos[0, 0] == 5.0
