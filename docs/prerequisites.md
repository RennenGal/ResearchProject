# Prerequisites

This document outlines the system requirements and setup instructions for the Protein Data Collector project.

## Quick Installation

For a quick automated setup, use our installation scripts:

### Linux/macOS
```bash
# Clone the repository
git clone https://github.com/RennenGal/ResearchProject.git
cd ResearchProject

# Run the installation script
./scripts/install.sh
```

### Windows (Command Prompt)
```cmd
REM Clone the repository
git clone https://github.com/RennenGal/ResearchProject.git
cd ResearchProject

REM Run the installation script
scripts\install.bat
```

### Windows (PowerShell)
```powershell
# Clone the repository
git clone https://github.com/RennenGal/ResearchProject.git
cd ResearchProject

# Run the installation script
.\scripts\install.ps1
```

**Note**: If you encounter PowerShell execution policy issues, run:
```powershell
Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser
```

## Manual Installation

If you prefer manual installation or the scripts don't work for your system, follow the detailed instructions below.

## System Requirements

### Minimum Requirements
- **Python**: 3.8 or higher (3.9+ recommended)
- **Memory**: 4GB RAM minimum, 8GB recommended
- **Storage**: 2GB free disk space
- **Network**: Internet connection for API access to UniProt and InterPro

### Supported Operating Systems
- **Linux**: Ubuntu 18.04+, CentOS 7+, or equivalent
- **macOS**: 10.15 (Catalina) or later
- **Windows**: Windows 10 or Windows 11

## Installation Instructions

### Linux (Ubuntu/Debian)

#### 1. Install Python 3.8+
```bash
# Update package list
sudo apt update

# Install Python 3.8+ and pip
sudo apt install python3.8 python3.8-venv python3.8-dev python3-pip

# Verify installation
python3.8 --version
pip3 --version
```

#### 2. Install Git
```bash
sudo apt install git
```

#### 3. Install Build Tools (for some Python packages)
```bash
sudo apt install build-essential libssl-dev libffi-dev
```

#### 4. Clone and Setup Project
```bash
# Clone the repository
git clone https://github.com/RennenGal/ResearchProject.git
cd ResearchProject

# Create virtual environment
python3.8 -m venv venv

# Activate virtual environment
source venv/bin/activate

# Upgrade pip
pip install --upgrade pip

# Install dependencies
pip install -r requirements.txt

# Install development dependencies (optional)
pip install -r requirements-dev.txt
```

### macOS

#### 1. Install Homebrew (if not already installed)
```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

#### 2. Install Python 3.8+
```bash
# Install Python via Homebrew
brew install python@3.9

# Verify installation
python3.9 --version
pip3 --version
```

#### 3. Install Git (if not already installed)
```bash
brew install git
```

#### 4. Clone and Setup Project
```bash
# Clone the repository
git clone https://github.com/RennenGal/ResearchProject.git
cd ResearchProject

# Create virtual environment
python3.9 -m venv venv

# Activate virtual environment
source venv/bin/activate

# Upgrade pip
pip install --upgrade pip

# Install dependencies
pip install -r requirements.txt

# Install development dependencies (optional)
pip install -r requirements-dev.txt
```

### Windows

#### 1. Install Python 3.8+
1. Download Python from [python.org](https://www.python.org/downloads/)
2. Run the installer and **check "Add Python to PATH"**
3. Verify installation in Command Prompt:
```cmd
python --version
pip --version
```

#### 2. Install Git
1. Download Git from [git-scm.com](https://git-scm.com/download/win)
2. Run the installer with default settings

#### 3. Clone and Setup Project
```cmd
# Clone the repository
git clone https://github.com/RennenGal/ResearchProject.git
cd ResearchProject

# Create virtual environment
python -m venv venv

# Activate virtual environment
venv\Scripts\activate

# Upgrade pip
python -m pip install --upgrade pip

# Install dependencies
pip install -r requirements.txt

# Install development dependencies (optional)
pip install -r requirements-dev.txt
```

## Database Setup

The project uses SQLite by default, which requires no additional setup. Database files are stored in the `db/` directory.

### Optional: MySQL Setup
If you prefer to use MySQL instead of SQLite:

#### Linux/macOS
```bash
# Install MySQL server
# Ubuntu/Debian:
sudo apt install mysql-server

# macOS:
brew install mysql
brew services start mysql
```

#### Windows
1. Download MySQL from [mysql.com](https://dev.mysql.com/downloads/mysql/)
2. Run the installer and follow the setup wizard

## Configuration

### 1. Environment Variables
Copy the example environment file and configure it:
```bash
cp .env.example .env
```

Edit `.env` with your preferred settings:
```bash
# API Configuration
UNIPROT_BASE_URL=https://rest.uniprot.org/
INTERPRO_BASE_URL=https://www.ebi.ac.uk/interpro/api/

# Rate Limiting
UNIPROT_REQUESTS_PER_SECOND=5.0
INTERPRO_REQUESTS_PER_SECOND=10.0

# Database (SQLite - default)
DATABASE_URL=sqlite:///db/protein_data.db

# Logging
LOG_LEVEL=INFO
```

### 2. Configuration Files
The project includes configuration files in the `config/` directory:
- `config/development.json` - Development settings
- `config/production.json` - Production settings

## Verification

### 1. Run Tests
```bash
# Activate virtual environment (if not already active)
source venv/bin/activate  # Linux/macOS
# or
venv\Scripts\activate     # Windows

# Run all tests
python -m pytest tests/ -v

# Expected output: All tests should pass
```

### 2. Test CLI Tool
```bash
# Test the command-line interface
python -m protein_data_collector.cli --help

# Test API connectivity
python -m protein_data_collector.cli test-apis
```

### 3. Verify Installation
```bash
# Check if the package is properly installed
python -c "import protein_data_collector; print('✅ Installation successful!')"
```

## Development Tools (Optional)

### Code Formatting and Linting
```bash
# Install pre-commit hooks (recommended for contributors)
pre-commit install

# Format code with Black
black protein_data_collector/ tests/

# Sort imports with isort
isort protein_data_collector/ tests/

# Run linting with flake8
flake8 protein_data_collector/ tests/

# Type checking with mypy
mypy protein_data_collector/
```

### IDE Setup
The project works well with:
- **VS Code**: Install Python extension
- **PyCharm**: Open project directory
- **Vim/Neovim**: Use with Python LSP

## Troubleshooting

### Common Issues

#### Python Version Issues
```bash
# Check Python version
python --version

# If using wrong version, specify explicitly
python3.9 -m venv venv  # Use specific version
```

#### Permission Issues (Linux/macOS)
```bash
# If pip install fails with permissions
pip install --user -r requirements.txt

# Or use virtual environment (recommended)
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

#### Windows Path Issues
- Ensure Python is added to PATH during installation
- Use `python` instead of `python3` on Windows
- Use `Scripts\activate` instead of `bin/activate`

#### Network Issues
```bash
# Test API connectivity
curl -I https://rest.uniprot.org/
curl -I https://www.ebi.ac.uk/interpro/api/

# If behind corporate firewall, configure proxy:
pip install --proxy http://proxy.company.com:port -r requirements.txt
```

### Getting Help
- Check the [README.md](../README.md) for project overview
- Review [installation.md](installation.md) for detailed setup
- Open an issue on GitHub for bugs or questions
- Check the `docs/` directory for additional documentation

## Next Steps

After completing the prerequisites:
1. Read the [README.md](../README.md) for project overview
2. Check [installation.md](installation.md) for detailed setup instructions
3. Review [scripts-workflow-guide.md](scripts-workflow-guide.md) for usage examples
4. Explore the `scripts/` directory for data collection tools

## Hardware Recommendations

### For Development
- **CPU**: 2+ cores
- **RAM**: 8GB
- **Storage**: SSD recommended for better performance

### For Production Data Collection
- **CPU**: 4+ cores (for concurrent API requests)
- **RAM**: 16GB+ (for large datasets)
- **Storage**: 50GB+ SSD (for database and logs)
- **Network**: Stable internet connection with good bandwidth

## Security Considerations

- Keep API keys secure (use environment variables)
- Regularly update dependencies: `pip install --upgrade -r requirements.txt`
- Use virtual environments to isolate dependencies
- Follow rate limiting guidelines for external APIs