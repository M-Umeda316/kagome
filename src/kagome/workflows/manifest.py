"""Run manifest for experiment tracking."""
from __future__ import annotations

import json
import subprocess
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


@dataclass
class RunManifest:
    """実験のプロベナンスを記録するマニフェスト。

    CLAUDE.md 非交渉要件: seed, config_path, git_sha, backend, output_dir を
    すべての実験で記録する。save() で JSON にシリアライズ。
    """

    config_path: str
    seed: int
    backend: str
    output_dir: str
    git_sha: str = ''
    git_dirty: bool = False
    timestamp: str = ''
    extra: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.timestamp:
            self.timestamp = datetime.now(timezone.utc).isoformat()
        if not self.git_sha:
            # Auto-resolve provenance: record the dirty flag alongside the SHA so a
            # recorded commit can be trusted to match the executed code (RF17).
            self.git_sha = _get_git_sha()
            self.git_dirty = _get_git_dirty()

    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(asdict(self), indent=2), encoding='utf-8')

    @classmethod
    def append_resume(
        cls, path: Path, ckpt_step: int, ckpt_cycle: int,
    ) -> bool:
        """Append a resume record to an existing manifest without clobbering it (W1).

        The original ``git_sha``/``timestamp``/``seed`` identify the run that
        produced the bulk of the trajectory; a resume only continues it and must
        not overwrite that provenance. The resuming code's git state is recorded
        under ``extra['resume_history']`` so both segments stay auditable.

        Returns True if the record was appended, False if *path* does not exist
        (the caller then falls back to a fresh :meth:`save`).
        """
        if not path.exists():
            return False
        data = json.loads(path.read_text(encoding='utf-8'))
        extra = data.setdefault('extra', {})
        history = extra.setdefault('resume_history', [])
        history.append({
            'timestamp': datetime.now(timezone.utc).isoformat(),
            'git_sha': _get_git_sha(),
            'git_dirty': _get_git_dirty(),
            'ckpt_step': int(ckpt_step),
            'ckpt_cycle': int(ckpt_cycle),
        })
        path.write_text(json.dumps(data, indent=2), encoding='utf-8')
        return True

    @staticmethod
    def record_production_start_step(path: Path, step: int) -> None:
        """Record the post-equilibration step where the production (cycle) loop
        begins, into ``extra['production_start_step']`` (A4).

        Activation and equilibration run outside :meth:`PolymerizationWorkflow.run`,
        so ``state.step`` at the start of the cycle loop is the true production
        onset that k_p fitting must anchor t=0 to. Called only on fresh runs; a
        resume preserves the original run's recorded value (do not overwrite).
        """
        if not path.exists():
            return
        data = json.loads(path.read_text(encoding='utf-8'))
        data.setdefault('extra', {})['production_start_step'] = int(step)
        path.write_text(json.dumps(data, indent=2), encoding='utf-8')


def _normalize_value(v: Any) -> Any:
    """Convert numpy scalars to Python built-ins for JSON serialisation."""
    if isinstance(v, dict):
        return {k: _normalize_value(val) for k, val in v.items()}
    if isinstance(v, (list, tuple)):
        return [_normalize_value(item) for item in v]
    if hasattr(v, 'item'):
        return v.item()
    return v


def _get_git_sha() -> str:
    try:
        result = subprocess.run(
            ['git', 'rev-parse', 'HEAD'],
            capture_output=True, text=True, timeout=5,
        )
        return result.stdout.strip() if result.returncode == 0 else 'unknown'
    except Exception:
        return 'unknown'


def _get_git_dirty() -> bool:
    """True if the working tree has uncommitted changes (tracked or untracked).

    A clean SHA is only a faithful identifier of the executed code when the tree
    is clean; record this so a dirty run is auditable (RF17).
    """
    try:
        result = subprocess.run(
            ['git', 'status', '--porcelain'],
            capture_output=True, text=True, timeout=5,
        )
        if result.returncode != 0:
            return False
        return bool(result.stdout.strip())
    except Exception:
        return False
