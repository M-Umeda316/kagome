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
    config_path: str
    seed: int
    backend: str
    output_dir: str
    git_sha: str = ''
    timestamp: str = ''
    extra: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.timestamp:
            self.timestamp = datetime.now(timezone.utc).isoformat()
        if not self.git_sha:
            self.git_sha = _get_git_sha()

    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(asdict(self), indent=2), encoding='utf-8')


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
