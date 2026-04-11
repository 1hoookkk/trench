param(
  [switch]$DryRun = $false
)

$ErrorActionPreference = "Stop"

$SourceDir  = "C:\Users\hooki\Trench"
$DestDir    = "C:\Users\hooki\trench-juce\plugin\assets"

# --- Whitelist: only these cross the boundary ---
$ApprovedBodies = @(
  "speaker_knockerz.json",
  "aluminum_siding.json",
  "small_talk.json",
  "cul_de_sac.json"
)

# --- 1. Sonic tables ---
$tablesSrc = "$SourceDir\docs\sonic_tables\tables.json"
$tablesDst = "$DestDir\sonic_tables.json"
if (Test-Path $tablesSrc) {
  if ($DryRun) { Write-Host "[dry-run] would copy: $tablesSrc -> $tablesDst" }
  else {
    Copy-Item $tablesSrc -Destination $tablesDst -Force
    Write-Host "OK  sonic_tables.json" -ForegroundColor Green
  }
} else {
  Write-Host "MISSING  $tablesSrc" -ForegroundColor Red
  exit 1
}

# --- 2. Shipping bodies (whitelist only) ---
$cartDst = "$DestDir\cartridges"
$missing = @()

foreach ($body in $ApprovedBodies) {
  $src = "$SourceDir\cartridges\$body"
  if (Test-Path $src) {
    if ($DryRun) { Write-Host "[dry-run] would copy: $src -> $cartDst\$body" }
    else {
      Copy-Item $src -Destination "$cartDst\$body" -Force
      Write-Host "OK  $body" -ForegroundColor Green
    }
  } else {
    Write-Host "MISSING  $body (not found at $src)" -ForegroundColor Red
    $missing += $body
  }
}

# --- 3. Warn about non-shipping bodies in dest ---
$existing = Get-ChildItem "$cartDst\*.json" | ForEach-Object { $_.Name }
$stowaways = $existing | Where-Object { $_ -notin $ApprovedBodies }

if ($stowaways.Count -gt 0) {
  Write-Host ""
  Write-Host "STOWAWAYS in plugin cartridges (not on whitelist):" -ForegroundColor Yellow
  foreach ($s in $stowaways) {
    Write-Host "  $s" -ForegroundColor Yellow
  }
  Write-Host "Remove these before tagging a release." -ForegroundColor Yellow
}

# --- Summary ---
Write-Host ""
if ($missing.Count -gt 0) {
  Write-Host "BLOCKED: $($missing.Count) shipping bodies not found at source." -ForegroundColor Red
  exit 1
} else {
  Write-Host "Sync complete." -ForegroundColor Green
}
