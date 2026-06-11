"""Tests for JSONL trajectory writer and reader."""
import json

import pytest

from src.io.trajectory import TrajectoryFrame, TrajectoryWriter
from src.io.readers import read_trajectory


def _make_frame(step: int = 0, cycle: int = 0, phase: str = 'biased') -> TrajectoryFrame:
    return TrajectoryFrame(
        step=step, time_fs=step * 0.25, phase=phase, cycle=cycle,
        energy_base=-10.0, energy_bias=5.0, energy_total=-5.0,
        positions=[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
    )


class TestTrajectoryWriter:

    def test_creates_file_with_header(self, tmp_path):
        path = tmp_path / 'traj.jsonl'
        writer = TrajectoryWriter(path, species=['C', 'C'], save_interval=1)
        writer.close()

        lines = path.read_text(encoding='utf-8').strip().split('\n')
        header = json.loads(lines[0])
        assert header['_header'] is True
        assert header['species'] == ['C', 'C']
        assert header['n_atoms'] == 2

    def test_writes_frames(self, tmp_path):
        path = tmp_path / 'traj.jsonl'
        writer = TrajectoryWriter(path, species=['C', 'C'], save_interval=1)
        writer.write_frame(_make_frame(step=0))
        writer.write_frame(_make_frame(step=1))
        writer.close()

        lines = path.read_text(encoding='utf-8').strip().split('\n')
        assert len(lines) == 3  # header + 2 frames
        assert writer.n_frames == 2

    def test_should_write_interval(self, tmp_path):
        path = tmp_path / 'traj.jsonl'
        writer = TrajectoryWriter(path, species=['C'], save_interval=10)

        assert writer.should_write(0) is True
        assert writer.should_write(5) is False
        assert writer.should_write(10) is True
        assert writer.should_write(20) is True
        writer.close()

    def test_should_write_disabled(self, tmp_path):
        path = tmp_path / 'traj.jsonl'
        writer = TrajectoryWriter(path, species=['C'], save_interval=0)
        assert writer.should_write(0) is False
        assert writer.should_write(100) is False
        writer.close()

    def test_metadata_in_header(self, tmp_path):
        path = tmp_path / 'traj.jsonl'
        writer = TrajectoryWriter(
            path, species=['C'], save_interval=1,
            metadata={'config': 'test', 'seed': 42},
        )
        writer.close()

        lines = path.read_text(encoding='utf-8').strip().split('\n')
        header = json.loads(lines[0])
        assert header['metadata']['seed'] == 42


class TestReadTrajectory:

    def test_roundtrip(self, tmp_path):
        path = tmp_path / 'traj.jsonl'
        writer = TrajectoryWriter(path, species=['C', 'O'], save_interval=1)
        f1 = _make_frame(step=0, cycle=0, phase='biased')
        f2 = _make_frame(step=100, cycle=0, phase='unbiased')
        writer.write_frame(f1)
        writer.write_frame(f2)
        writer.close()

        header, frames = read_trajectory(path)

        assert header['species'] == ['C', 'O']
        assert len(frames) == 2
        assert frames[0].step == 0
        assert frames[0].phase == 'biased'
        assert frames[1].step == 100
        assert frames[1].energy_total == -5.0

    def test_empty_trajectory(self, tmp_path):
        path = tmp_path / 'traj.jsonl'
        writer = TrajectoryWriter(path, species=['H'], save_interval=1)
        writer.close()

        header, frames = read_trajectory(path)
        assert header['n_atoms'] == 1
        assert len(frames) == 0
