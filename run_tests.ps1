$ErrorActionPreference = "Stop"
Set-Location -Path $PSScriptRoot
$env:PYTHONPATH = Join-Path $PSScriptRoot "src"
$vpy = Join-Path $PSScriptRoot ".venv\Scripts\python.exe"
if (-not (Test-Path $vpy)) { throw "Run .\setup.ps1 first." }
& $vpy -m pytest -q tests\test_orthogonal_validation_post5.py
Write-Host "--- full suite ---" -ForegroundColor Cyan
& $vpy -m pytest -q
