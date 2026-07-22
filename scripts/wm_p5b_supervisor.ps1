# WM-P5b campaign supervisor (Windows PowerShell 5.1)
#
# Problem this script solves: the WM-P5b GPU campaign (20 runs, driven by
# scripts/run_wm_p5b_campaign.sh inside WSL pfpoly-gpu) is periodically hit by
# a known INFRASTRUCTURE fault, not a code bug -- after roughly 11h of
# sustained orb/CUDA load, the WSL GPU/CUDA driver context degrades and every
# subsequent orb run fails with `torch.AcceleratorError: CUDA error: unknown
# error` (confirmed: a no-mixing run failed identically, and a tiny CUDA
# matmul still succeeds -- so it is the driver/VM state, not kagome or the
# mixing code). The only known remedy is `wsl --shutdown`, which tears down
# the whole WSL VM and resets the CUDA context; runs succeed again afterwards.
#
# This script is a Windows-side supervisor that drives the 20-run matrix one
# (seed, mix_ps) pair at a time via scripts/run_wm_p5b_campaign.sh, and
# PROACTIVELY resets WSL every -Batch successful runs so the campaign never
# has to ride out the ~11h degradation window uncontrolled. It never
# implements the run itself -- scripts/run_wm_p5b_campaign.sh remains the
# single source of truth for the matrix, flags, retry policy, and
# idempotency (SKIP-DONE via summary.json). See that script's header and
# specs/decisions.md "追補 2026-07-19 -- WM-P5b 実験計画 (sweep x 多シード,
# ユーザー承認スコープ C)" for the run design this supervisor wraps.
#
# Usage:
#   .\scripts\wm_p5b_supervisor.ps1                          # real campaign
#   .\scripts\wm_p5b_supervisor.ps1 -DryRun                  # exercise control flow only, no wsl/GPU touched
#   .\scripts\wm_p5b_supervisor.ps1 -DryRun -DryRunFailSet @('s7_mix25')
#   .\scripts\wm_p5b_supervisor.ps1 -Batch 3 -MaxAttemptsPerRun 4
#
# Safe to launch detached (e.g. Start-Process / a scheduled task) for the
# ~2.5 day campaign duration; every action is timestamped to both the
# console and runs\wm_p5b\supervisor.log so progress can be audited after
# the fact even if the interactive session is gone.
#
# --- DryRun design note (ambiguity resolution) -----------------------------
# The brief's single `-DryRunFailSet` example needs to express two distinct
# simulated behaviours across the validation scenarios: (a) a run that fails
# once and then recovers after a WSL reset, and (b) a run that keeps failing
# even after a reset (the persistent/escalation case). One flat list cannot
# encode both without ambiguity, so this script splits it into two
# independently injectable sets, keeping the brief's example name for the
# recoverable case since that is the common/expected transient-GPU-fault
# shape in production:
#   -DryRunFailSet       : fails on its FIRST simulated attempt only, then
#                          succeeds on the post-reset retry (models a
#                          transient CUDA fault that the reset actually
#                          fixes -- the normal case this supervisor exists
#                          for).
#   -DryRunAlwaysFailSet : fails on every simulated attempt, including after
#                          a reset (models a persistent/non-GPU problem that
#                          Reset-Wsl cannot fix -- exercises the
#                          two-consecutive-failures ABORT path).
# -DryRunPreDoneSet seeds Test-RunDone as already-true for the listed run
# names before the loop starts, simulating runs that already have a valid
# summary.json on disk (idempotent restart scenario).

[CmdletBinding()]
param(
    [int]$Batch = 2,
    [int]$MaxAttemptsPerRun = 3,
    [string]$Distro = 'Ubuntu-24.04',
    [string]$RepoWsl = '/mnt/c/Users/shanu/Documents/Python/kagome',
    [switch]$DryRun,
    [string[]]$DryRunFailSet = @(),
    [string[]]$DryRunAlwaysFailSet = @(),
    [string[]]$DryRunPreDoneSet = @()
)

$ErrorActionPreference = 'Stop'

# Force WSL to emit UTF-8 without the UTF-16LE console mangling that can
# otherwise garble `wsl -d ... -- echo ...` output under Windows PowerShell
# 5.1's default console codepage handling.
$env:WSL_UTF8 = '1'

$RepoRoot = Split-Path -Parent $PSScriptRoot
$RunsDir  = Join-Path $RepoRoot 'runs\wm_p5b'
New-Item -ItemType Directory -Force -Path $RunsDir | Out-Null
# NOTE: was 'supervisor.log' — that file got wedged with a stale Windows
# write-lock (a WSL `tail -f` on a Windows-drive log left a dangling handle
# that survived wsl --shutdown), so logging is redirected to a fresh file.
$LogPath  = Join-Path $RunsDir 'supervisor-2.log'

# --- DryRun simulated state (only ever touched when -DryRun is set) --------
$script:DryRunState    = @{}   # run name -> $true (done) / $false (not done)
$script:DryRunAttempts = @{}   # run name -> number of Invoke-OneRun calls so far
$script:ResetCount     = 0

foreach ($n in $DryRunPreDoneSet) { $script:DryRunState[$n] = $true }

function Write-Log {
    param([Parameter(Mandatory)][string]$Message)
    $ts   = Get-Date -Format 'yyyy-MM-dd HH:mm:ss'
    $line = "[$ts] $Message"
    Write-Host $line
    Add-Content -Path $LogPath -Value $line
}

function Test-RunDone {
    param([Parameter(Mandatory)][int]$Seed, [Parameter(Mandatory)][int]$Mix)
    $name = 's{0}_mix{1}' -f $Seed, $Mix

    if ($DryRun) {
        if ($script:DryRunState.ContainsKey($name)) { return $script:DryRunState[$name] }
        return $false
    }

    $summaryPath = Join-Path $RunsDir "$name\summary.json"
    if (-not (Test-Path -LiteralPath $summaryPath)) { return $false }
    $item = Get-Item -LiteralPath $summaryPath -ErrorAction SilentlyContinue
    if (-not $item -or $item.Length -eq 0) { return $false }
    try {
        $null = Get-Content -Raw -LiteralPath $summaryPath | ConvertFrom-Json
        return $true
    } catch {
        return $false
    }
}

function Invoke-WslProbe {
    # Best-effort `wsl -d $Distro -- <args>` invocation for probes. Never
    # throws -- a transient probe error just yields empty output so the
    # caller's retry loop treats it as "not ready yet". Real wsl calls are
    # skipped entirely in -DryRun.
    #
    # TIMEOUT (added after a wedge): right after `wsl --shutdown`, a `wsl`
    # probe can hang FOREVER, not on wsl.exe (which exits) but on
    # `... | Out-String` blocking because a lingering WSL-side process still
    # holds the pipe's write end. That froze the whole supervisor for 40+ min
    # during a proactive reset. We therefore run the probe in a background job
    # and abandon it after $TimeoutSec, returning '' (== "not ready yet") so
    # the caller's retry loop advances instead of blocking indefinitely.
    param([Parameter(Mandatory)][string[]]$Args, [int]$TimeoutSec = 30)
    $job = Start-Job -ScriptBlock {
        param($d, $a)
        & wsl -d $d -- @a 2>&1 | Out-String
    } -ArgumentList $Distro, $Args
    if (Wait-Job $job -Timeout $TimeoutSec) {
        $out = Receive-Job $job
        Remove-Job $job -Force -ErrorAction SilentlyContinue
        return ("$out").Trim()
    }
    # Timed out: stop waiting on the (possibly wedged) probe and move on.
    Stop-Job $job -ErrorAction SilentlyContinue
    Remove-Job $job -Force -ErrorAction SilentlyContinue
    Write-Log "  (probe timed out after ${TimeoutSec}s, treating as not-ready)"
    return ''
}

function Reset-Wsl {
    param([string]$Reason = 'scheduled')
    $script:ResetCount++
    Write-Log "RESET-WSL #$($script:ResetCount) ($Reason): issuing 'wsl --shutdown' to clear the degraded GPU/CUDA context..."

    if ($DryRun) {
        Write-Log '  [DryRun] would run: wsl --shutdown'
    } else {
        try {
            & wsl --shutdown 2>&1 | Out-Null
        } catch {
            Write-Log "  WARNING: 'wsl --shutdown' raised an error, continuing to probe anyway: $($_.Exception.Message)"
        }
    }

    Write-Log '  Waiting 15s for the WSL VM to fully tear down...'
    Start-Sleep -Seconds 15

    Write-Log "  Probing for WSL ($Distro) to come back up..."
    $wslUp = $false
    for ($i = 1; $i -le 10; $i++) {
        if ($DryRun) {
            $out = 'wsl-up'
        } else {
            $out = Invoke-WslProbe -Args @('echo', 'wsl-up')
        }
        if ($out -match 'wsl-up') {
            $wslUp = $true
            Write-Log "  WSL is back up (probe $i/10)."
            break
        }
        Write-Log "  WSL not yet responsive (probe $i/10): '$out'"
        Start-Sleep -Seconds 5
    }
    if (-not $wslUp) {
        throw 'WSL did not return after shutdown -- manual intervention required'
    }

    Write-Log '  Waiting for GPU memory to drain (<2000 MiB)...'
    $drained = $false
    for ($i = 1; $i -le 24; $i++) {
        if ($DryRun) {
            $used = 0
        } else {
            $raw = Invoke-WslProbe -Args @('nvidia-smi', '--query-gpu=memory.used', '--format=csv,noheader,nounits')
            $firstLine = ($raw -split "`n" | Select-Object -First 1).Trim()
            $parsed = 0
            if ([int]::TryParse($firstLine, [ref]$parsed)) { $used = $parsed } else { $used = 999999 }
        }
        if ($used -lt 2000) {
            $drained = $true
            Write-Log "  GPU drained ($used MiB used) after $i probe(s)."
            break
        }
        Start-Sleep -Seconds 5
    }
    if (-not $drained) {
        Write-Log '  WARNING: GPU did not report <2000 MiB within the probe budget; proceeding anyway (non-fatal).'
    }
}

function Invoke-OneRun {
    # Runs exactly one (seed, mix) pair via scripts/run_wm_p5b_campaign.sh,
    # blocking until it returns. Success is judged by the CALLER via
    # Test-RunDone (summary.json), never by this function's exit code --
    # the launcher can FAIL-STUCK/FAIL-EXHAUSTED internally and still exit
    # in a way that isn't a reliable success/failure signal on its own.
    param([Parameter(Mandatory)][int]$Seed, [Parameter(Mandatory)][int]$Mix)
    $name = 's{0}_mix{1}' -f $Seed, $Mix
    $cmdStr = "cd $RepoWsl && SEEDS='$Seed' MIX_PS='$Mix' MAX_ATTEMPTS=$MaxAttemptsPerRun bash scripts/run_wm_p5b_campaign.sh >> runs/wm_p5b/campaign.out 2>&1"

    Write-Log "INVOKE $name : wsl -d $Distro -- bash -lc `"$cmdStr`""

    if ($DryRun) {
        if (-not $script:DryRunAttempts.ContainsKey($name)) { $script:DryRunAttempts[$name] = 0 }
        $script:DryRunAttempts[$name]++
        $attemptNum = $script:DryRunAttempts[$name]

        if ($DryRunAlwaysFailSet -contains $name) {
            $script:DryRunState[$name] = $false
            Write-Log "  [DryRun] simulated outcome for $name : FAIL (persistent, attempt $attemptNum)"
        } elseif (($DryRunFailSet -contains $name) -and ($attemptNum -eq 1)) {
            $script:DryRunState[$name] = $false
            Write-Log "  [DryRun] simulated outcome for $name : FAIL (transient, attempt $attemptNum)"
        } else {
            $script:DryRunState[$name] = $true
            Write-Log "  [DryRun] simulated outcome for $name : OK (attempt $attemptNum)"
        }
        return
    }

    try {
        # $cmdStr is passed as a single pre-built argument, not re-quoted
        # inline -- this sidesteps PowerShell 5.1's backtick/quote escaping
        # entirely, since PowerShell hands the string to wsl.exe as one
        # argv entry and bash -lc receives it verbatim.
        & wsl -d $Distro -- bash -lc $cmdStr
    } catch {
        Write-Log "  WARNING: wsl invocation for $name raised an error (ignored -- judged via summary.json, not exit code): $($_.Exception.Message)"
    }
    Write-Log "  Launcher call for $name returned (LASTEXITCODE=$LASTEXITCODE, informational only)."
}

# --- Build the fixed 20-run matrix: mix_ps outer, seed inner ---------------
# (identical enumeration order to scripts/run_wm_p5b_campaign.sh's default
# MIX_PS/SEEDS loops, so supervisor and launcher logs line up run-for-run.)
$MixList  = @(0, 25, 50, 100)
$SeedList = @(7, 11, 17, 23, 42)
$Pairs = foreach ($mix in $MixList) {
    foreach ($seed in $SeedList) {
        [PSCustomObject]@{ Seed = $seed; Mix = $mix; Name = 's{0}_mix{1}' -f $seed, $mix }
    }
}

Write-Log '=== WM-P5b supervisor start ==='
Write-Log "  DryRun:            $($DryRun.IsPresent)"
Write-Log "  Batch:             $Batch"
Write-Log "  MaxAttemptsPerRun: $MaxAttemptsPerRun"
Write-Log "  Distro:            $Distro"
Write-Log "  RepoWsl:           $RepoWsl"
if ($DryRun) {
    Write-Log "  DryRunFailSet:       $($DryRunFailSet -join ', ')"
    Write-Log "  DryRunAlwaysFailSet: $($DryRunAlwaysFailSet -join ', ')"
    Write-Log "  DryRunPreDoneSet:    $($DryRunPreDoneSet -join ', ')"
}

# The GPU is known to be currently degraded going into this campaign, so we
# always start with a fresh reset regardless of how many runs are already done.
Reset-Wsl -Reason 'initial (GPU assumed degraded at supervisor start)'
$sinceReset      = 0
$consecutiveFail = 0
$aborted         = $false

foreach ($pair in $Pairs) {
    $seed = $pair.Seed
    $mix  = $pair.Mix
    $name = $pair.Name

    if (Test-RunDone -Seed $seed -Mix $mix) {
        Write-Log "SKIP-DONE $name (valid summary.json already present)"
        continue
    }

    if ($sinceReset -ge $Batch) {
        Reset-Wsl -Reason "proactive (after $sinceReset successful runs, Batch=$Batch)"
        $sinceReset = 0
    }

    Invoke-OneRun -Seed $seed -Mix $mix

    if (Test-RunDone -Seed $seed -Mix $mix) {
        Write-Log "OK $name"
        $sinceReset++
        $consecutiveFail = 0
        continue
    }

    Write-Log "RUN-FAILED $name (likely GPU degraded) -> reset + one retry"
    Reset-Wsl -Reason "reactive (after $name failed)"
    $sinceReset = 0
    Invoke-OneRun -Seed $seed -Mix $mix

    if (Test-RunDone -Seed $seed -Mix $mix) {
        Write-Log "OK-after-reset $name"
        # Counts toward the next proactive-reset batch just like a
        # first-attempt OK -- otherwise the Batch=N cadence would silently
        # drift after every recovered run, undermining the whole point of
        # proactive resets (see script header note on this ambiguity).
        $sinceReset++
        $consecutiveFail = 0
    } else {
        Write-Log "FAILED-AFTER-RESET $name"
        $consecutiveFail++
        if ($consecutiveFail -ge 2) {
            Write-Log 'ABORT: 2 consecutive runs failed even after a fresh WSL reset -- persistent problem, escalating'
            $aborted = $true
            break
        }
    }
}

$doneNames   = @($Pairs | Where-Object { Test-RunDone -Seed $_.Seed -Mix $_.Mix } | ForEach-Object { $_.Name })
$notDoneNames = @($Pairs | Where-Object { -not (Test-RunDone -Seed $_.Seed -Mix $_.Mix) } | ForEach-Object { $_.Name })

Write-Log '=== WM-P5b supervisor summary ==='
Write-Log "  Done:       $($doneNames.Count)/$($Pairs.Count)"
Write-Log "  Resets:     $($script:ResetCount)"
Write-Log "  Aborted:    $aborted"
if ($notDoneNames.Count -gt 0) {
    Write-Log "  Not done:   $($notDoneNames -join ', ')"
}
Write-Log '=== WM-P5b supervisor end ==='
