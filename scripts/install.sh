#!/bin/bash
# Installation script for Linux/macOS - Protein Data Collector
# This script automates the setup process for the Protein Data Collector project

set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Function to check if command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to detect OS
detect_os() {
    if [[ "$OSTYPE" == "linux-gnu"* ]]; then
        echo "linux"
    elif [[ "$OSTYPE" == "darwin"* ]]; then
        echo "macos"
    else
        echo "unknown"
    fi
}

# Function to check Python version
check_python_version() {
    local python_cmd=$1
    if command_exists "$python_cmd"; then
        local version=$($python_cmd --version 2>&1 | cut -d' ' -f2)
        local major=$(echo $version | cut -d'.' -f1)
        local minor=$(echo $version | cut -d'.' -f2)
        
        if [ "$major" -eq 3 ] && [ "$minor" -ge 8 ]; then
            echo "$python_cmd"
            return 0
        fi
    fi
    return 1
}

# Function to find suitable Python
find_python() {
    local python_commands=("python3.11" "python3.10" "python3.9" "python3.8" "python3" "python")
    
    for cmd in "${python_commands[@]}"; do
        if check_python_version "$cmd"; then
            return 0
        fi
    done
    return 1
}

# Main installation function
main() {
    print_status "Starting Protein Data Collector installation..."
    
    # Detect OS
    OS=$(detect_os)
    print_status "Detected OS: $OS"
    
    # Check if we're in the right directory
    if [ ! -f "pyproject.toml" ] || [ ! -f "requirements.txt" ]; then
        print_error "Please run this script from the project root directory"
        exit 1
    fi
    
    # Find Python
    print_status "Checking Python installation..."
    if ! PYTHON_CMD=$(find_python); then
        print_error "Python 3.8+ is required but not found"
        print_status "Please install Python 3.8+ and try again"
        
        if [ "$OS" = "linux" ]; then
            print_status "On Ubuntu/Debian: sudo apt install python3.9 python3.9-venv python3.9-dev"
        elif [ "$OS" = "macos" ]; then
            print_status "On macOS: brew install python@3.9"
        fi
        exit 1
    fi
    
    print_success "Found Python: $PYTHON_CMD ($($PYTHON_CMD --version))"
    
    # Check Git
    print_status "Checking Git installation..."
    if ! command_exists git; then
        print_error "Git is required but not found"
        if [ "$OS" = "linux" ]; then
            print_status "Install with: sudo apt install git"
        elif [ "$OS" = "macos" ]; then
            print_status "Install with: brew install git"
        fi
        exit 1
    fi
    print_success "Git is installed"
    
    # Install system dependencies (Linux only)
    if [ "$OS" = "linux" ]; then
        print_status "Installing system dependencies..."
        if command_exists apt-get; then
            sudo apt-get update
            sudo apt-get install -y build-essential libssl-dev libffi-dev python3-dev
        elif command_exists yum; then
            sudo yum groupinstall -y "Development Tools"
            sudo yum install -y openssl-devel libffi-devel python3-devel
        else
            print_warning "Could not detect package manager. You may need to install build tools manually."
        fi
    fi
    
    # Create virtual environment
    print_status "Creating virtual environment..."
    if [ -d "venv" ]; then
        print_warning "Virtual environment already exists. Removing old one..."
        rm -rf venv
    fi
    
    $PYTHON_CMD -m venv venv
    print_success "Virtual environment created"
    
    # Activate virtual environment
    print_status "Activating virtual environment..."
    source venv/bin/activate
    
    # Upgrade pip
    print_status "Upgrading pip..."
    python -m pip install --upgrade pip
    
    # Install dependencies
    print_status "Installing Python dependencies..."
    pip install -r requirements.txt
    
    # Install development dependencies (optional)
    if [ -f "requirements-dev.txt" ]; then
        read -p "Install development dependencies? (y/N): " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            print_status "Installing development dependencies..."
            pip install -r requirements-dev.txt
            print_success "Development dependencies installed"
        fi
    fi
    
    # Create .env file if it doesn't exist
    if [ ! -f ".env" ]; then
        print_status "Creating .env file from template..."
        cp .env.example .env
        print_success ".env file created"
        print_warning "Please review and update .env file with your settings"
    fi
    
    # Create database directory
    print_status "Setting up database directory..."
    mkdir -p db
    
    # Run tests to verify installation
    print_status "Running tests to verify installation..."
    if python -m pytest tests/ -q --tb=no; then
        print_success "All tests passed!"
    else
        print_warning "Some tests failed. Installation may still work, but please check the output."
    fi
    
    # Test CLI
    print_status "Testing CLI tool..."
    if python -m protein_data_collector.cli --help > /dev/null 2>&1; then
        print_success "CLI tool is working"
    else
        print_warning "CLI tool test failed"
    fi
    
    # Installation complete
    echo
    print_success "🎉 Installation completed successfully!"
    echo
    print_status "Next steps:"
    echo "  1. Activate the virtual environment: source venv/bin/activate"
    echo "  2. Review and update .env file with your settings"
    echo "  3. Run tests: python -m pytest tests/"
    echo "  4. Test CLI: python -m protein_data_collector.cli --help"
    echo "  5. Check documentation in docs/ directory"
    echo
    print_status "For more information, see docs/prerequisites.md"
}

# Run main function
main "$@"