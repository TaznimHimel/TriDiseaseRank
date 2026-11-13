@echo off
setlocal
cd /d "%~dp0"
where Rscript >nul 2>nul
if errorlevel 1 (
  echo Rscript was not found. Install R first, then run setup.bat.
  pause
  exit /b 1
)
echo Starting TriDiseaseRank-X...
Rscript scripts\check_environment.R
if errorlevel 1 (
  echo.
  echo Environment check failed. Run setup.bat and review the error.
  pause
  exit /b 1
)
Rscript scripts\launch_app.R
if errorlevel 1 (
  echo.
  echo The app stopped with an error. Review the messages above.
  pause
  exit /b 1
)
pause
