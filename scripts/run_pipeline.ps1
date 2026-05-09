# One-shot: create venv if missing, install deps, run analysis.
# Usage from project root (spaceflight-damage-response):
#   powershell -ExecutionPolicy Bypass -File .\scripts\run_pipeline.ps1

$ErrorActionPreference = "Stop"
$Root = Split-Path -Parent $PSScriptRoot
Set-Location $Root

Write-Host "Working directory: $Root"

$venvPy = Join-Path $Root ".venv\Scripts\python.exe"
if (-not (Test-Path $venvPy)) {
    Write-Host "Creating virtual environment..."
    python -m venv .venv
}

& (Join-Path $Root ".venv\Scripts\Activate.ps1")

Write-Host "Installing dependencies..."
pip install -r requirements.txt -q

Write-Host "Running main.py..."
python main.py
