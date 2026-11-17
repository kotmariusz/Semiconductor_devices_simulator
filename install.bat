@echo off
rem SemiSim 3D - Windows installer: double-click this file.
rem Finds Python 3.9+ (64-bit) and runs install_win.py, which makes a private
rem environment (.venv), gets the solver core, runs a self-test and adds
rem Start-menu and Desktop shortcuts.  No administrator rights needed.
setlocal
cd /d "%~dp0"
title SemiSim 3D - install
set "PYEXE="
rem the Python launcher (python.org) first, then python on PATH; the Microsoft
rem Store stub of python.exe fails this test, so it is not used
py -3 -c "import sys; sys.exit(0 if sys.version_info >= (3, 9) else 1)" >nul 2>&1 && set "PYEXE=py -3"
if not defined PYEXE (
    python -c "import sys; sys.exit(0 if sys.version_info >= (3, 9) else 1)" >nul 2>&1 && set "PYEXE=python"
)
if not defined PYEXE (
    echo.
    echo Python 3.9 or newer ^(64-bit^) is needed and was not found.
    echo Get it from https://www.python.org/downloads/windows/ - in its installer
    echo tick "Add python.exe to PATH" - then double-click install.bat again.
    echo.
    pause
    exit /b 1
)
%PYEXE% "%~dp0install_win.py" %*
set "RC=%ERRORLEVEL%"
echo.
if "%RC%"=="0" (
    echo SemiSim 3D is installed: start it from the Start menu or the Desktop shortcut.
) else (
    echo The installation did not finish - see the messages above.
)
pause
exit /b %RC%
