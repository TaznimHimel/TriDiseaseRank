@echo off
setlocal
cd /d "%~dp0"
where Rscript >nul 2>nul
if errorlevel 1 (
  echo Rscript was not found.
  echo Install R for Windows 64-bit first: https://cran.r-project.org/bin/windows/base/
  pause
  exit /b 1
)
echo Setting up TriDiseaseRank-X. This may take a while...
Rscript scripts\setup_windows.R
if errorlevel 1 (
  echo.
  echo Setup did not complete. Read the messages above and setup_package_report.csv.
  pause
  exit /b 1
)
echo.
echo Setup complete. You can now double-click run_app.bat.
pause
