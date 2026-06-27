"""Runtime diagnostics for long MD loops (no effect on physics or results).

Added while investigating the ``runs/s6_full_200x10`` cycle-15 hang: a long
paper-scale run that stalled with *normal* energy and temperature and no NaN —
the signature of a backend/GPU stall (most likely VRAM exhaustion -> WSL
system-memory fallback), not a TDBB/physics bug. See the
``PYTORCH_CUDA_ALLOC_CONF`` note in ``scripts/run_vinyl_aibn.py`` and
``specs/decisions.md`` (2026-06-15 VRAM record, 2026-06-27 cycle-15 hang).

:class:`StepWatchdog` wraps each MD step (arm() before, step_done() after) to:

* **re-arm faulthandler** so a step that hangs longer than ``KAGOME_WATCHDOG_S``
  dumps *every* thread's stack — this pinpoints whether the stall is inside
  ``calculator.compute`` / warp-CUDA or in Python logic, which a post-mortem of
  a killed process cannot;
* **WARN on a slow step** (> ``KAGOME_STALL_WARN_S``) with that step's peak
  VRAM, so a transient memory/cost spike is caught even if it recedes;
* **heartbeat** every ``KAGOME_HEARTBEAT_STEPS`` steps (step rate + VRAM) to
  confirm liveness and show how far a run got before stalling.

The thresholds are diagnostic knobs (not scientific hyperparameters) and default
to values that stay silent in normal operation, so any run is self-diagnosing
without code changes. The watchdog uses a *single-shot* timer re-armed each step
rather than a periodic dump, so it fires only on a genuine stall.
"""
from __future__ import annotations

import faulthandler
import logging
import os
import time

logger = logging.getLogger(__name__)


class StepWatchdog:
    """Per-step stall watchdog + VRAM/heartbeat logger. Diagnostic only.

    Does not touch positions, velocities, forces, RNG or any simulation state.
    """

    def __init__(self) -> None:
        self.watchdog_s = float(os.environ.get('KAGOME_WATCHDOG_S', '180'))
        self.warn_s = float(os.environ.get('KAGOME_STALL_WARN_S', '20'))
        self.heartbeat_steps = int(os.environ.get('KAGOME_HEARTBEAT_STEPS', '0'))
        self._tic = 0.0
        self._win_tic: float | None = None
        self._count = 0
        # Detect a CUDA device once; VRAM logging is skipped on CPU/MACE runs.
        self._torch = None
        self._cuda = False
        try:
            import torch
            self._torch = torch
            self._cuda = bool(torch.cuda.is_available())
        except Exception:
            pass
        if self.watchdog_s > 0 and not faulthandler.is_enabled():
            # Also surfaces a traceback on a fatal signal (e.g. SIGSEGV) — cheap.
            faulthandler.enable()

    def _vram_note(self) -> str:
        if not self._cuda:
            return ''
        gib = 1024 ** 3
        # torch-only view: max_memory_allocated captures THIS step's transient
        # backward-graph spike (peak reset in arm()). But warp/orb allocate
        # OUTSIDE torch's caching allocator, so these undercount the true device
        # footprint. mem_get_info() reads free/total straight from the CUDA
        # driver and so includes warp's memory — this is what actually governs
        # the WSL system-memory fallback that presents as a hang.
        peak = self._torch.cuda.max_memory_allocated() / gib
        note = f' (torch_peak={peak:.2f} GiB'
        try:
            free_b, total_b = self._torch.cuda.mem_get_info()
            used = (total_b - free_b) / gib
            note += f'; device used={used:.2f}/{total_b / gib:.2f} free={free_b / gib:.2f} GiB'
        except Exception:
            pass
        return note + ')'

    def arm(self) -> None:
        """Call immediately before an MD step."""
        if self.watchdog_s > 0:
            faulthandler.dump_traceback_later(self.watchdog_s, repeat=False)
        if self._cuda:
            # Reset so the peak read in step_done() belongs to THIS step, making a
            # transient spike attributable even after it recedes.
            self._torch.cuda.reset_peak_memory_stats()
        self._tic = time.perf_counter()
        if self._win_tic is None:
            self._win_tic = self._tic

    def step_done(self, *, phase: str, cycle: int, step: int) -> None:
        """Call immediately after an MD step completes."""
        now = time.perf_counter()
        if self.watchdog_s > 0:
            faulthandler.cancel_dump_traceback_later()
        dt = now - self._tic
        if self.warn_s > 0 and dt > self.warn_s:
            logger.warning(
                'Slow MD step: %.1f s at %s cycle %d step %d%s',
                dt, phase, cycle, step, self._vram_note(),
            )
        self._count += 1
        if self.heartbeat_steps > 0 and self._count % self.heartbeat_steps == 0:
            elapsed = now - (self._win_tic or now)
            rate = self.heartbeat_steps / elapsed if elapsed > 0 else float('nan')
            logger.info(
                'Heartbeat: %s cycle %d step %d, %.2f steps/s%s',
                phase, cycle, step, rate, self._vram_note(),
            )
            self._win_tic = now
