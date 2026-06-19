"""Tests for JSONL trajectory writer and reader."""
import json

import pytest

from src.io.trajectory import TrajectoryFrame, TrajectoryWriter
from src.io.readers import read_trajectory


def _make_frame(step: int = 0, cycle: int = 0, phase: str = 'biased',
                temperature_K: float = 0.0) -> TrajectoryFrame:
    return TrajectoryFrame(
        step=step, time_fs=step * 0.25, phase=phase, cycle=cycle,
        energy_base=-10.0, energy_bias=5.0, energy_total=-5.0,
        positions=[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        temperature_K=temperature_K,
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

    def test_temperature_roundtrip(self, tmp_path):
        path = tmp_path / 'traj.jsonl'
        writer = TrajectoryWriter(path, species=['C', 'C'], save_interval=1)
        writer.write_frame(_make_frame(step=0, temperature_K=312.5))
        writer.write_frame(_make_frame(step=1, temperature_K=487.3))
        writer.close()

        _, frames = read_trajectory(path)
        assert frames[0].temperature_K == pytest.approx(312.5)
        assert frames[1].temperature_K == pytest.approx(487.3)

    def test_temperature_defaults_to_zero(self, tmp_path):
        """Old trajectories without temperature_K field are read without error."""
        import json
        path = tmp_path / 'traj.jsonl'
        header = {'_header': True, 'species': ['C'], 'n_atoms': 1, 'save_interval': 1}
        frame = {
            'step': 0, 'time_fs': 0.0, 'phase': 'biased', 'cycle': 0,
            'energy_base': 0.0, 'energy_bias': 0.0, 'energy_total': 0.0,
            'positions': [[0.0, 0.0, 0.0]],
        }
        path.write_text(
            json.dumps(header) + '\n' + json.dumps(frame) + '\n',
            encoding='utf-8',
        )
        _, frames = read_trajectory(path)
        assert frames[0].temperature_K == 0.0

    def test_schema_version_in_header(self, tmp_path):
        path = tmp_path / 'traj.jsonl'
        writer = TrajectoryWriter(path, species=['C', 'O'], save_interval=1)
        writer.close()

        header, _ = read_trajectory(path)
        assert header['schema_version'] == 1

    def test_old_header_without_schema_version(self, tmp_path):
        """Pre-v1 trajectories without schema_version are readable (version 0)."""
        import json
        path = tmp_path / 'traj.jsonl'
        header = {'_header': True, 'species': ['C'], 'n_atoms': 1, 'save_interval': 1}
        path.write_text(json.dumps(header) + '\n', encoding='utf-8')

        loaded_header, _ = read_trajectory(path)
        assert loaded_header.get('schema_version', 0) == 0
