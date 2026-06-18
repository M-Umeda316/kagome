# S6: Paper-scale vinyl radical polymerization reproduction (Windows PowerShell)
# Paper: arXiv:2511.22874, Table S1, SI S-3
# System: 200 methyl acrylate monomers + 10 AIBN initiators (approx 2520 atoms)
# Required GPU VRAM: >=24 GB (OrbMol-v2 single-step footprint ~9.5 GB at this scale)
#
# Usage (from repo root, with pfpoly-gpu conda env active):
#   .\scripts\run_s6_paper_scale.ps1
#   .\scripts\run_s6_paper_scale.ps1 -Seed 42 -OutputDir runs\s6_seed42
#   .\scripts\run_s6_paper_scale.ps1 -NCycles 100
#
# For 16 GB GPU (half-scale 100+5), run run_vinyl_aibn.py directly:
#   python scripts\run_vinyl_aibn.py --seed 7 --output-dir runs\s6_half_scale `
#       --n-monomers 100 --n-initiators 5 --activation --activation-f2 0.3 `
#       --activation-f1-max 250 --f2 5.0 --density 0.5 --temperature 333.0 `
#       --no-barostat --backend orb --device cuda --n-cycles 30 `
#       --biased-steps 2000 --unbiased-steps 500 --equil-steps 2000 --timestep-fs 1.0
#
# --- Parameter rationale ---
# n_monomers=200, n_initiators=10  : Paper Table S1, Section 3 (200+10 system)
# activation                        : AIBN V^d C-N dissociation before propagation
# activation-f2=0.3, f1_max=250    : OrbMol-v2 C-N barrier ~39 kcal/mol
#                                     See specs/decisions.md 2026-06-18
# f2=5.0                           : Validated capture radius for OrbMol-v2 PES
#                                     See specs/decisions.md 2026-06-17
# n_cycles=50                      : Paper reports multi-hundred cycles; 50 is a feasible start
# biased_steps=2000                : Validated in S2-S3 runs
# unbiased_steps=500               : Validated in S2-S3 runs
# timestep_fs=1.0                  : Standard for organic ML MD
# density=0.5                      : Paper SI S-3 (methyl acrylate melt)
# temperature=333.0                : Paper Table S1 (60 C)
# no-barostat                      : NVT validated
# seed=7                           : Matches S2-S3 best-result seed

param(
    [int]$Seed = 7,
    [string]$OutputDir = "",
    [string]$Device = "cuda",
    [int]$NCycles = 50,
    [int]$BiasedSteps = 2000,
    [int]$UnbiasedSteps = 500,
    [int]$EquilSteps = 2000
)

$ErrorActionPreference = "Stop"

if ($OutputDir -eq "") {
    $OutputDir = "runs\s6_paper_scale_seed$Seed"
}

# Ensure VRAM fragmentation workaround is set
if (-not $env:PYTORCH_CUDA_ALLOC_CONF) {
    $env:PYTORCH_CUDA_ALLOC_CONF = "expandable_segments:True"
}
if (-not $env:KMP_DUPLICATE_LIB_OK) {
    $env:KMP_DUPLICATE_LIB_OK = "TRUE"
}

Write-Host "=== S6 paper-scale run ==="
Write-Host "  Seed:           $Seed"
Write-Host "  Output dir:     $OutputDir"
Write-Host "  Device:         $Device"
Write-Host "  N cycles:       $NCycles"
Write-Host "  Biased steps:   $BiasedSteps"
Write-Host "  Unbiased steps: $UnbiasedSteps"
Write-Host ""

# VRAM check
try {
    $nvidiaSmi = nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits 2>$null
    if ($nvidiaSmi) {
        $vramMB = [int]($nvidiaSmi.Trim().Split("`n")[0])
        $vramGB = [math]::Round($vramMB / 1024, 1)
        Write-Host "  GPU VRAM:       $vramGB GB (need >=24 GB for 200+10 system)"
        if ($vramMB -lt 20000) {
            Write-Host "  WARNING: VRAM may be insufficient. Consider --n-monomers 100 --n-initiators 5 for 16 GB GPU."
        }
    }
} catch {
    Write-Host "  GPU VRAM:       (nvidia-smi not available)"
}

Write-Host ""
Write-Host "Starting run..."

python scripts/run_vinyl_aibn.py `
    --seed $Seed `
    --output-dir $OutputDir `
    --n-monomers 200 `
    --n-initiators 10 `
    --activation `
    --activation-f2 0.3 `
    --activation-f1-max 250 `
    --activation-steps 5000 `
    --f2 5.0 `
    --density 0.5 `
    --temperature 333.0 `
    --no-barostat `
    --backend orb `
    --device $Device `
    --n-cycles $NCycles `
    --biased-steps $BiasedSteps `
    --unbiased-steps $UnbiasedSteps `
    --equil-steps $EquilSteps `
    --timestep-fs 1.0 `
    --minimize `
    --minimize-fmax 1.0

if ($LASTEXITCODE -ne 0) {
    Write-Host "ERROR: run_vinyl_aibn.py failed with exit code $LASTEXITCODE"
    exit $LASTEXITCODE
}

Write-Host ""
Write-Host "Run complete. Generating figures..."

# n_reactive_sites is auto-read from the trajectory header (RF2: alpha denominator
# = n_monomers, written by PolymerizationWorkflow). No CLI override needed.
python scripts/reproduce_figures.py `
    --trajectory "$OutputDir\trajectory.jsonl" `
    --bonds "$OutputDir\bonds.jsonl" `
    --target-temperature 333.0 `
    --timestep-fs 1.0 `
    --output-dir "$OutputDir\figures"

Write-Host ""
Write-Host "=== S6 done ==="
Write-Host "  Artifacts: $OutputDir\"
Write-Host "  Figures:   $OutputDir\figures\"
Write-Host ""
Write-Host "Reproduction command:"
Write-Host "  .\scripts\run_s6_paper_scale.ps1 -Seed $Seed -OutputDir $OutputDir"
