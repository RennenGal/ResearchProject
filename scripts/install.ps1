# Installation script for Windows PowerShell - Protein Data Collector
# This script automates the setup process for the Protein Data Collector project

# Set error action preference
$ErrorActionPreference = "Stop"

# Function to write colored output
function Write-Status {
    param([string]$Message)
    Write-Host "[INFO] $Message" -ForegroundColor Blue
}

function Write-Success {
    param([string]$Message)
    Write-Host "[SUCCESS] $Message" -ForegroundColor Green
}

function Write-Warning {
    param([string]$Message)
    Write-Host "[WARNING] $Message" -ForegroundColor Yellow
}

function Write-Error {
    param([string]$Message)
    Write-Host "[ERROR] $Message" -ForegroundColor Red
}

# Function to check if command exists
function Test-Command {
    param([string]$Command)
    try {
        Get-Command $Command -ErrorAction Stop | Out-Null
        return $true
    }
    catch {
        return $false
    }
}

# Function to check Python version
function Test-PythonVersion {
    param([string]$PythonCommand)
    
    try {
        $version = & $PythonCommand --version 2>&1
        if ($version -match "Python (\d+)\.(\d+)\.(\d+)") {
            $major = [int]$matches[1]
            $minor = [int]$matches[2]
            
            if ($major -eq 3 -and $minor -ge 8) {
                return $true
            }
        }
        return $false
    }
    catch {
        return $false
    }
}

# Main installation function
function Install-ProteinDataCollector {
    Write-Status "Starting Protein Data Collector installation..."
    
    # Check if we're in the right directory
    if (-not (Test-Path "pyproject.toml") -or -not (Test-Path "requirements.txt")) {
        Write-Error "Please run this script from the project root directory"
        Read-Host "Press Enter to exit"
        exit 1
    }
    
    # Check Python installation
    Write-Status "Checking Python installation..."
    $pythonCommands = @("python", "python3", "py")
    $pythonCmd = $null
    
    foreach ($cmd in $pythonCommands) {
        if (Test-Command $cmd) {
            if (Test-PythonVersion $cmd) {
                $pythonCmd = $cmd
                break
            }
        }
    }
    
    if (-not $pythonCmd) {
        Write-Error "Python 3.8+ is required but not found"
        Write-Status "Please install Python 3.8+ from https://python.org/downloads/"
        Write-Status "Make sure to check 'Add Python to PATH' during installation"
        Read-Host "Press Enter to exit"
        exit 1
    }
    
    $pythonVersion = & $pythonCmd --version
    Write-Success "Found Python: $pythonVersion"
    
    # Check Git installation
    Write-Status "Checking Git installation..."
    if (-not (Test-Command "git")) {
        Write-Error "Git is required but not found"
        Write-Status "Please install Git from https://git-scm.com/download/win"
        Read-Host "Press Enter to exit"
        exit 1
    }
    Write-Success "Git is installed"
    
    # Remove existing virtual environment if it exists
    if (Test-Path "venv") {
        Write-Warning "Virtual environment already exists. Removing old one..."
        Remove-Item -Recurse -Force "venv"
    }
    
    # Create virtual environment
    Write-Status "Creating virtual environment..."
    try {
        & $pythonCmd -m venv venv
        Write-Success "Virtual environment created"
    }
    catch {
        Write-Error "Failed to create virtual environment: $_"
        Read-Host "Press Enter to exit"
        exit 1
    }
    
    # Activate virtual environment
    Write-Status "Activating virtual environment..."
    $activateScript = "venv\Scripts\Activate.ps1"
    
    # Check if we can run the activation script
    try {
        & $activateScript
    }
    catch {
        Write-Warning "Could not activate virtual environment with PowerShell script"
        Write-Status "You may need to run: Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser"
        Write-Status "Trying alternative activation method..."
        
        # Use the batch file instead
        cmd /c "venv\Scripts\activate.bat && python -m pip install --upgrade pip"
    }
    
    # Upgrade pip
    Write-Status "Upgrading pip..."
    try {
        python -m pip install --upgrade pip
    }
    catch {
        Write-Warning "Failed to upgrade pip, continuing anyway..."
    }
    
    # Install dependencies
    Write-Status "Installing Python dependencies..."
    try {
        pip install -r requirements.txt
        Write-Success "Dependencies installed"
    }
    catch {
        Write-Error "Failed to install dependencies: $_"
        Read-Host "Press Enter to exit"
        exit 1
    }
    
    # Install development dependencies (optional)
    if (Test-Path "requirements-dev.txt") {
        $installDev = Read-Host "Install development dependencies? (y/N)"
        if ($installDev -eq "y" -or $installDev -eq "Y") {
            Write-Status "Installing development dependencies..."
            try {
                pip install -r requirements-dev.txt
                Write-Success "Development dependencies installed"
            }
            catch {
                Write-Warning "Failed to install some development dependencies: $_"
            }
        }
    }
    
    # Create .env file if it doesn't exist
    if (-not (Test-Path ".env")) {
        Write-Status "Creating .env file from template..."
        try {
            Copy-Item ".env.example" ".env"
            Write-Success ".env file created"
            Write-Warning "Please review and update .env file with your settings"
        }
        catch {
            Write-Warning "Failed to create .env file: $_"
        }
    }
    
    # Create database directory
    Write-Status "Setting up database directory..."
    if (-not (Test-Path "db")) {
        New-Item -ItemType Directory -Path "db" | Out-Null
    }
    
    # Run tests to verify installation
    Write-Status "Running tests to verify installation..."
    try {
        python -m pytest tests/ -q --tb=no | Out-Null
        Write-Success "All tests passed!"
    }
    catch {
        Write-Warning "Some tests failed. Installation may still work, but please check manually."
    }
    
    # Test CLI
    Write-Status "Testing CLI tool..."
    try {
        python -m protein_data_collector.cli --help | Out-Null
        Write-Success "CLI tool is working"
    }
    catch {
        Write-Warning "CLI tool test failed"
    }
    
    # Installation complete
    Write-Host ""
    Write-Success "🎉 Installation completed successfully!"
    Write-Host ""
    Write-Status "Next steps:"
    Write-Host "  1. Activate the virtual environment: venv\Scripts\Activate.ps1"
    Write-Host "  2. Review and update .env file with your settings"
    Write-Host "  3. Run tests: python -m pytest tests/"
    Write-Host "  4. Test CLI: python -m protein_data_collector.cli --help"
    Write-Host "  5. Check documentation in docs\ directory"
    Write-Host ""
    Write-Status "For more information, see docs\prerequisites.md"
    Write-Host ""
    Read-Host "Press Enter to exit"
}

# Run the installation
try {
    Install-ProteinDataCollector
}
catch {
    Write-Error "Installation failed: $_"
    Read-Host "Press Enter to exit"
    exit 1
}