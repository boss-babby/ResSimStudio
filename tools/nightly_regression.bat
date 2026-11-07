@echo off
rem C33 nightly deck regression (Windows). Register once with:
rem   schtasks /Create /SC DAILY /ST 02:00 /TN "ResSim nightly regression" /TR "\"%~f0\""
rem Needs Python 3 (for the download step) and ResSimCLI.exe next to this file or in RESSIM_HOME.
setlocal
if "%RESSIM_HOME%"=="" set RESSIM_HOME=%~dp0
if "%RESSIM_TESTS%"=="" set RESSIM_TESTS=%RESSIM_HOME%opm-tests
cd /d "%RESSIM_HOME%"
python "%RESSIM_HOME%tools\fetch_opm_tests.py" -o "%RESSIM_TESTS%" --update
for /f "tokens=1-3 delims=/-. " %%a in ("%date%") do set STAMP=%%c%%b%%a
"%RESSIM_HOME%ResSimCLI.exe" --regress "%RESSIM_TESTS%\manifest.txt" -o "%RESSIM_HOME%regression\%STAMP%" --history "%RESSIM_HOME%regression\history.tsv" > "%RESSIM_HOME%regression\last_run.log" 2>&1
if errorlevel 1 (echo REGRESSIONS FOUND - see regression\%STAMP%\report.html) else (echo all cases pass)
endlocal
