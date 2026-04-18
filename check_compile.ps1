$cl = 'C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.44.35207\bin\Hostx64\x64\cl.exe'
$src = 'C:\Users\hooki\trench-juce\plugin\source\TrenchEngine.cpp'
$includes = @(
  'C:\Users\hooki\trench-juce\plugin\build\TRENCH_artefacts\JuceLibraryCode',
  'C:\Users\hooki\trench-juce\plugin\juce\modules',
  'C:\Users\hooki\trench-juce\plugin\build\JUCE_BinaryData_Assets\JuceLibraryCode',
  'C:\Users\hooki\trench-juce\plugin\modules',
  'C:\Users\hooki\trench-juce\plugin\source'
)
$iflags = ($includes | ForEach-Object { "/I`"$_`"" }) -join ' '
$defs = '/D _MBCS /D WIN32 /D _WINDOWS /D NDEBUG /D JUCE_SHARED_CODE=1 /D JUCE_MODULE_AVAILABLE_juce_core=1 /D JUCE_GLOBAL_MODULE_SETTINGS_INCLUDED=1'
$allargs = "/c /nologo /W4 /WX- /EHsc /MT /fp:fast /Zc:wchar_t /Zc:forScope /Zc:inline /std:c++latest /Zc:__cplusplus /TP $defs $iflags /Fo`"C:\Users\hooki\Trench\trenchengine_check.obj`" `"$src`""
Write-Host "Compiling TrenchEngine.cpp..."
$psi = New-Object System.Diagnostics.ProcessStartInfo
$psi.FileName = $cl
$psi.Arguments = $allargs
$psi.UseShellExecute = $false
$psi.RedirectStandardOutput = $true
$psi.RedirectStandardError = $true
$proc = [System.Diagnostics.Process]::Start($psi)
$stdout = $proc.StandardOutput.ReadToEnd()
$stderr = $proc.StandardError.ReadToEnd()
$proc.WaitForExit()
Write-Host $stdout
Write-Host $stderr
Write-Host "Exit code: $($proc.ExitCode)"
