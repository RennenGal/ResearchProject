@echo off
REM Installation script for Windows - Protein Data Collector
REM This script automates the setup process for the Protein Data Collector project

setlocal enabledelayedexpansion

REM Colors for output (limited in batch)
set "INFO=[INFO]"
set "SUCCESS=[SUCCESS]"
set "WARNING=[WARNING]"
set "ERROR=[ERROR]"

echo %INFO% Starting Protein Data Collector installation...

REM Check if we're in the right directory
if not exist "pyproject.toml" (
    echo %ERROR% Please run this script from the project root directory
    pause
    exit /b 1
)

if not exist "requirements.txt" (
    echo %ERROR% requirements.txt not found. Please run from project root directory
    pause
    exit /b 1
)

REM Check Python installation
echo %INFO% Checking Python installation...
python --version >nul 2>&1
if errorlevel 1 (
    echo %ERROR% Python is not installed or not in PATH
    echo %INFO% Please install Python 3.8+ from https://python.org/downloads/
    echo %INFO% Make sure to check "Add Python to PATH" during installation
    pause
    exit /b 1
)

REM Check Python version
for /f "tokens=2" %%i in ('python --version 2^>^&1') do set PYTHON_VERSION=%%i
echo %SUCCESS% Found Python: %PYTHON_VERSION%

REM Extract major and minor version numbers
for /f "tokens=1,2 delims=." %%a in ("%PYTHON_VERSION%") do (
    set MAJOR=%%a
    set MINOR=%%b
)

REM Check if version is 3.8+
if %MAJOR% LSS 3 (
    echo %ERROR% Python 3.8+ is required, found %PYTHON_VERSION%
    pause
    exit /b 1
)

if %MAJOR% EQU 3 if %MINOR% LSS 8 (
    echo %ERROR% Python 3.8+ is required, found %PYTHON_VERSION%
    pause
    exit /b 1
)

REM Check Git installation
echo %INFO% Checking Git installation...
git --version >nul 2>&1
if errorlevel 1 (
    echo %ERROR% Git is not installed or not in PATH
    echo %INFO% Please install Git from https://git-scm.com/download/win
    pause
    exit /b 1
)
echo %SUCCESS% Git is installed

REM Remove existing virtual environment if it exists
if exist "venv" (
    echo %WARNING% Virtual environment already exists. Removing old one...
    rmdir /s /q venv
)

REM Create virtual environment
echo %INFO% Creating virtual environment...
python -m venv venv
if errorlevel 1 (
    echo %ERROR% Failed to create virtual environment
    pause
    exit /b 1
)
echo %SUCCESS% Virtual environment created

REM Activate virtual environment
echo %INFO% Activating virtual environment...
call venv\Scripts\activate.bat

REM Upgrade pip
echo %INFO% Upgrading pip...
python -m pip install --upgrade pip
if errorlevel 1 (
    echo %WARNING% Failed to upgrade pip, continuing anyway...
)

REM Install dependencies
echo %INFO% Installing Python dependencies...
pip install -r requirements.txt
if errorlevel 1 (
    echo %ERROR% Failed to install dependencies
    pause
    exit /b 1
)
echo %SUCCESS% Dependencies installed

REM Install development dependencies (optional)
if exist "requirements-dev.txt" (
    set /p INSTALL_DEV="Install development dependencies? (y/N): "
    if /i "!INSTALL_DEV!"=="y" (
        echo %INFO% Installing development dependencies...
        pip install -r requirements-dev.txt
        if errorlevel 1 (
            echo %WARNING% Failed to install some development dependencies
        ) else (
            echo %SUCCESS% Development dependencies installed
        )
    )
)

REM Create .env file if it doesn't exist
if not exist ".env" (
    echo %INFO% Creating .env file from template...
    copy .env.example .env >nul
    if errorlevel 1 (
        echo %WARNING% Failed to create .env file
    ) else (
        echo %SUCCESS% .env file created
        echo %WARNING% Please review and update .env file with your settings
    )
)

REM Create database directory
echo %INFO% Setting up database directory...
if not exist "db" mkdir db

REM Run tests to verify installation
echo %INFO% Running tests to verify installation...
python -m pytest tests/ -q --tb=no >nul 2>&1
if errorlevel 1 (
    echo %WARNING% Some tests failed. Installation may still work, but please check manually.
) else (
    echo %SUCCESS% All tests passed!
)

REM Test CLI
echo %INFO% Testing CLI tool...
python -m protein_data_collector.cli --help >nul 2>&1
if errorlevel 1 (
    echo %WARNING% CLI tool test failed
) else (
    echo %SUCCESS% CLI tool is working
)

REM Installation complete
echo.
echo %SUCCESS% 🎉 Installation completed successfully!
echo.
echo %INFO% Next steps:
echo   1. Activate the virtual environment: venv\Scripts\activate
echo   2. Review and update .env file with your settings
echo   3. Run tests: python -m pytest tests/
echo   4. Test CLI: python -m protein_data_collector.cli --help
echo   5. Check documentation in docs\ directory
echo.
echo %INFO% For more information, see docs\prerequisites.md
echo.
pause