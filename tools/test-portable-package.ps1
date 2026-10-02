param([Parameter(Mandatory=$true)][string]$Archive)
$ErrorActionPreference = 'Stop'
# A fresh extraction also catches accidental dependence on the build directory.
$archivePath = (Resolve-Path -LiteralPath $Archive).Path
$testRoot = Join-Path (Split-Path $archivePath) ('smoke-' + [guid]::NewGuid().ToString('N'))
Expand-Archive -LiteralPath $archivePath -DestinationPath $testRoot
$executables = @(Get-ChildItem -LiteralPath $testRoot -Recurse -Filter NavierGui.exe)
if ($executables.Count -ne 1) { throw 'Expected exactly one packaged NavierGui.exe' }
$exe = $executables[0].FullName
$packageRoot = Split-Path (Split-Path $exe)
foreach ($relative in @('bin/Qt6Core.dll', 'bin/Qt6Gui.dll', 'bin/Qt6Widgets.dll',
    'bin/libgcc_s_seh-1.dll', 'bin/libstdc++-6.dll', 'bin/libwinpthread-1.dll',
    'bin/qt.conf', 'plugins/platforms/qwindows.dll', 'README.md')) {
    if (-not (Test-Path -LiteralPath (Join-Path $packageRoot $relative) -PathType Leaf)) {
        throw "Package is missing $relative"
    }
}
foreach ($relative in @('include', 'lib', 'tests')) {
    if (Test-Path -LiteralPath (Join-Path $packageRoot $relative)) {
        throw "Development files unexpectedly packaged: $relative"
    }
}
if (Get-ChildItem -LiteralPath $packageRoot -Recurse -File | Where-Object { $_.Extension -in @('.db', '.sqlite', '.sqlite3') }) {
    throw 'An experiment database was unexpectedly packaged'
}
$savedEnvironment = @{}
$variables = @('PATH', 'QT_PLUGIN_PATH', 'QT_QPA_PLATFORM_PLUGIN_PATH', 'QT_QPA_PLATFORM', 'QML_IMPORT_PATH', 'QML2_IMPORT_PATH')
try {
    foreach ($name in $variables) {
        $savedEnvironment[$name] = [Environment]::GetEnvironmentVariable($name, 'Process')
        [Environment]::SetEnvironmentVariable($name, $null, 'Process')
    }
    $env:PATH = "$env:SystemRoot\System32;$env:SystemRoot"
    # Exercise the real Windows plugin, not the offscreen test plugin.
    $env:QT_QPA_PLATFORM = 'windows'
    $process = Start-Process -FilePath $exe -ArgumentList '--smoke-test' -WorkingDirectory $testRoot -WindowStyle Hidden -PassThru
    if (-not $process.WaitForExit(15000)) {
        $process.Kill()
        throw 'Packaged GUI did not exit within 15 seconds'
    }
    if ($process.ExitCode -ne 0) { throw "Packaged GUI exited with $($process.ExitCode)" }
    Write-Output "Packaged Windows startup passed with Qt/compiler paths removed: $exe"
} finally {
    foreach ($name in $variables) {
        [Environment]::SetEnvironmentVariable($name, $savedEnvironment[$name], 'Process')
    }
}
# Retain the unique extraction for inspection; no existing files are deleted.
