param(
  [ValidateSet("shipping", "p2k_refs")]
  [string]$Set = "shipping",

  [string]$PackDir = "vault\\_vocal_pack_2026-04-11"
)

$ErrorActionPreference = "Stop"

$repoRoot = Resolve-Path (Join-Path $PSScriptRoot "..")
$livePath = Join-Path $repoRoot "trench_live.json"
$packRoot = Join-Path $repoRoot $PackDir
$setDir = Join-Path $packRoot $Set

if (-not (Test-Path $setDir)) {
  throw "Pack set folder not found: $setDir"
}

$files = Get-ChildItem -Path $setDir -Filter "*.json" | Sort-Object Name
if (-not $files -or $files.Count -lt 1) {
  throw "No JSON bodies found in: $setDir"
}

Write-Host ""
Write-Host "TRENCH · Hotload vocal pack"
Write-Host "Repo : $repoRoot"
Write-Host "Pack : $setDir"
Write-Host "Live : $livePath"
Write-Host ""

for ($i = 0; $i -lt $files.Count; $i++) {
  $n = $i + 1
  Write-Host ("{0,2}. {1}" -f $n, $files[$i].Name)
}

Write-Host ""
$pick = Read-Host "Pick a number (1-$($files.Count))"
$parsed = 0
if (-not [int]::TryParse($pick, [ref]$parsed)) {
  throw "Invalid selection: '$pick'"
}
$idx = $parsed - 1
if ($idx -lt 0 -or $idx -ge $files.Count) {
  throw "Selection out of range: '$pick'"
}

$src = $files[$idx].FullName
Copy-Item -Path $src -Destination $livePath -Force

Write-Host ""
Write-Host "Loaded:"
Write-Host "  $src"
Write-Host "→ $livePath"
Write-Host ""
