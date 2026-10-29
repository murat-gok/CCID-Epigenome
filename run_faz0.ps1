# FAZ 0 diagnostic gate
$ErrorActionPreference = "Stop"
Set-Location -Path $PSScriptRoot
$env:PYTHONPATH = Join-Path $PSScriptRoot "src"
$vpy = Join-Path $PSScriptRoot ".venv\Scripts\python.exe"
if (-not (Test-Path $vpy)) { throw "Run .\setup.ps1 first." }

$stem = "results\real\safe_osk_target_signature_O1O2"
$allLoci = "${stem}_all_loci.csv.gz"
$donorEff = "${stem}_donor_effects.csv.gz"
$sig = "${stem}.csv.gz"

$missing = @()
foreach ($f in @($allLoci, $donorEff)) { if (-not (Test-Path $f)) { $missing += $f } }
if ($missing.Count -gt 0) {
    Write-Host "Required discovery outputs are missing:" -ForegroundColor Red
    $missing | ForEach-Object { Write-Host "  $_" }
    Write-Host "`nThese are side outputs of scripts\03_build_target_signature.py." -ForegroundColor Yellow
    Write-Host "Follow START_HERE.md section 3, then rerun this script." -ForegroundColor Yellow
    Write-Host "`nIf your discovery output uses a different stem, locate it with:" -ForegroundColor DarkGray
    Write-Host '    Get-ChildItem -Recurse results -Filter *_all_loci.csv.gz'
    exit 1
}

$cliArgs = @(
    "scripts\faz0_diagnostics.py",
    "--all-loci", $allLoci,
    "--donor-effects", $donorEff,
    "--shortlist", "data\reference\frozen_primary_shortlist_v090.csv",
    "--outdir", "results\faz0"
)
if (Test-Path $sig) { $cliArgs += @("--signature", $sig) }
else { Write-Host "[warn] signature file not found; the frozen-list sanity check will be skipped." -ForegroundColor Yellow }

& $vpy @cliArgs
Write-Host "`nSend these four back for the route decision:" -ForegroundColor Cyan
Write-Host "    results\faz0\faz0_summary.json"
Write-Host "    results\faz0\testA_panel_floor_anatomy.csv"
Write-Host "    results\faz0\testB_ranking_sensitivity.csv"
Write-Host "    results\faz0\testC_global_delta_distribution.csv"
