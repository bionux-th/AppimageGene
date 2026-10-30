# Run on a real Windows desktop after building; keeps user settings untouched.
$ErrorActionPreference = 'Stop'
$Binary = Join-Path $PSScriptRoot '..\target\x86_64-pc-windows-msvc\release\appimagegene.exe'
if (-not (Test-Path $Binary)) { throw "Build first: cargo build --release -p hx-app --target x86_64-pc-windows-msvc" }
$OriginalAppData = $env:APPDATA
$Private = Join-Path ([IO.Path]::GetTempPath()) ('ag-smoke-' + [Guid]::NewGuid().ToString())
New-Item -ItemType Directory -Path (Join-Path $Private 'appimagegene') | Out-Null
try {
    $env:APPDATA = $Private
    [IO.File]::WriteAllText((Join-Path $Private 'appimagegene\preferences'),"language=zh-CN`n",(New-Object Text.UTF8Encoding $false))
    Write-Host "Private configuration: $Private"
    Write-Host 'Check native Open/Save/Folder dialogs, Chinese labels, Explorer reveal and PDF fallback. Close AG when finished.'
    Start-Process -FilePath $Binary -Wait
    Get-ChildItem (Join-Path $Private 'appimagegene')
} finally {
    $env:APPDATA = $OriginalAppData
    Write-Host "Smoke-test files retained at $Private (no original configuration changed)."
}
