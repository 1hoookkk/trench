param(
  [string]$HostAddr = "127.0.0.1",
  [int]$Port = 8000,
  [switch]$Reload = $false,
  [switch]$NoBrowser = $false
)

$ErrorActionPreference = "Stop"

$repoRoot = Resolve-Path (Join-Path $PSScriptRoot "..")
Set-Location $repoRoot

$baseUrl = "http://$HostAddr`:$Port"

Write-Host ""
Write-Host "TRENCH authoring runtime"
Write-Host "Repo: $repoRoot"
Write-Host ""
Write-Host "Surfaces:"
Write-Host "  Terrain  : $baseUrl/terrain"
Write-Host "  Designer : $baseUrl/"
Write-Host "  Author   : $baseUrl/author"
Write-Host "  Workbench: $baseUrl/workbench"
Write-Host "  Phonemes : $baseUrl/phonemes"
Write-Host "  VaultLab : $baseUrl/vaultlab"
Write-Host ""
Write-Host "Tip: Terrain writes to ``trench_live.json`` when the gate passes."
Write-Host ""

if (-not $NoBrowser) {
  Start-Process "$baseUrl/terrain" | Out-Null
}

$uvicornArgs = @(
  "-m", "uvicorn", "pyruntime.api:app",
  "--host", $HostAddr,
  "--port", "$Port"
)

if ($Reload) {
  $uvicornArgs += @(
    "--reload",
    "--reload-dir", "pyruntime",
    "--reload-exclude", "vault/**",
    "--reload-exclude", "target/**",
    "--reload-exclude", "target-*/**",
    "--reload-exclude", ".git/**",
    "--reload-exclude", "datasets/**",
    "--reload-exclude", "cartridges/**",
    "--reload-exclude", ".codex/**",
    "--reload-exclude", ".claude/**",
    "--reload-exclude", ".vscode/**",
    "--reload-delay", "1.0"
  )
}

& python @uvicornArgs
