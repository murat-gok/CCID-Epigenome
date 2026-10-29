# CCID-Epigenome environment bootstrap (Windows PowerShell)
$ErrorActionPreference = "Stop"
Set-Location -Path $PSScriptRoot

Write-Host "=== CCID-Epigenome setup ===" -ForegroundColor Cyan

# locate a Python launcher
$pyExe = $null
$pyArgs = @()
if (Get-Command py -ErrorAction SilentlyContinue) {
    $pyExe = "py"; $pyArgs = @("-3")
} elseif (Get-Command python -ErrorAction SilentlyContinue) {
    $pyExe = "python"
} else {
    throw "Python 3.11+ was not found on PATH. Install it from python.org, tick 'Add python.exe to PATH', then reopen PowerShell."
}
Write-Host ("Python launcher: " + $pyExe + " " + ($pyArgs -join " "))
& $pyExe @pyArgs --version

if (-not (Test-Path ".venv")) {
    Write-Host "Creating .venv ..." -ForegroundColor Yellow
    & $pyExe @pyArgs -m venv .venv
}

$vpy = Join-Path $PSScriptRoot ".venv\Scripts\python.exe"
if (-not (Test-Path $vpy)) { throw "Virtual environment creation failed: $vpy not found" }

Write-Host "Installing dependencies ..." -ForegroundColor Yellow
& $vpy -m pip install --upgrade pip --quiet
& $vpy -m pip install -r requirements.txt --quiet
if ($LASTEXITCODE -ne 0) { throw "Dependency installation failed" }

Write-Host "`nRunning test suite ..." -ForegroundColor Yellow
$env:PYTHONPATH = Join-Path $PSScriptRoot "src"
& $vpy -m pytest -q
$code = $LASTEXITCODE

Write-Host ""
if ($code -eq 0) {
    Write-Host "Setup complete. Expect '72 passed'." -ForegroundColor Green
    Write-Host "Next: START_HERE.md section 3 — regenerate the FAZ 0 inputs." -ForegroundColor Green
} else {
    Write-Host "Not all tests passed. Send the output before continuing." -ForegroundColor Red
}

Write-Host "`nActivate the environment in new shells with:" -ForegroundColor Cyan
Write-Host "    .\.venv\Scripts\Activate.ps1"
Write-Host "If activation is blocked, run once:" -ForegroundColor DarkGray
Write-Host "    Set-ExecutionPolicy -Scope CurrentUser RemoteSigned"
