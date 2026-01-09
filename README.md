# Protein Data Collector

A comprehensive Python application for collecting and analyzing protein family data from InterPro, with a focus on TIM barrel protein structures.

## Overview

This project provides tools to collect, store, and analyze protein family data from the InterPro database. It specializes in identifying and cataloging TIM (Triosephosphate Isomerase) barrel protein structures, which are important enzyme folds found across many biological processes.

## Key Features

- **Unified Data Collection**: Collects both PFAM families and InterPro entries with TIM barrel annotations
- **Hybrid Search Strategy**: Uses comprehensive search terms to find all relevant TIM barrel entries
- **Database Integration**: Stores data in MySQL with proper schema and relationships
- **API Integration**: Interfaces with InterPro REST API with rate limiting and caching
- **Comprehensive Testing**: Includes unit tests and integration tests
- **Clean Architecture**: Modular design with separation of concerns

## Current Status

✅ **49 TIM Barrel Entries Collected**
- 18 PFAM families (PF entries)
- 31 InterPro entries (IPR entries)
- Unified storage in single database table
- No artificial limits - collects all available entries

## Quick Start

### Prerequisites

- Python 3.8+
- MySQL database
- Virtual environment (recommended)

### Installation

1. Clone the repository:
```bash
git clone https://github.com/RennenGal/ResearchProject.git
cd ResearchProject
```

2. Create and activate virtual environment:
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. Install dependencies:
```bash
pip install -r requirements.txt
```

4. Set up configuration:
```bash
cp config.json.example config.test.json
# Edit config.test.json with your database settings
```

5. Set up database:
```bash
# Create MySQL database and user
mysql -u root -p < test_schema.sql
```

### Usage

#### Collect TIM Barrel Data
```bash
# Production collection
python scripts/collect_tim_barrel_entries.py

# Dry run (test without storing)
python scripts/collect_tim_barrel_entries.py --dry-run

# Verbose output
python scripts/collect_tim_barrel_entries.py --verbose
```

#### View Collection Status
```bash
python scripts/tim_barrel_summary.py
```

#### Run Tests
```bash
pytest
```

## Project Structure

```
├── protein_data_collector/     # Main application package
│   ├── api/                   # API clients (InterPro, UniProt)
│   ├── database/              # Database models and connections
│   ├── models/                # Data models and validation
│   ├── collector/             # Data collection logic
│   └── query/                 # Query processing
├── scripts/                   # Collection and utility scripts
├── tests/                     # Test suite
├── config/                    # Configuration files
└── docs/                      # Documentation (if any)
```

## Architecture

### Data Collection Strategy

The application uses a **hybrid search approach** to ensure comprehensive coverage:

1. **Phase 1: Direct PFAM Family Search**
   - Searches PFAM families using multiple TIM barrel-related terms
   - Filters results to ensure relevance
   - Captures protein family classifications

2. **Phase 2: InterPro Entry Search**
   - Searches InterPro entries (IPR records) for structural classifications
   - Includes domains, superfamilies, and active sites
   - Captures broader structural annotations

3. **Unified Storage**
   - All entries stored in single `tim_barrel_entries` table
   - Automatic deduplication by accession
   - Preserves entry type and metadata

### Database Schema

- **tim_barrel_entries**: Unified table for both PFAM and InterPro entries
- **interpro_proteins**: Protein sequences with InterPro annotations
- Proper foreign key relationships and indexing

## Configuration

The application uses JSON configuration files:

- `config.test.json`: Test environment settings
- `config/development.json`: Development settings
- `config/production.json`: Production settings

Key configuration sections:
- Database connection settings
- API rate limiting
- Logging configuration
- Cache settings

## Testing

Comprehensive test suite includes:
- Unit tests for individual components
- Integration tests for API interactions
- Database tests with fixtures
- 92 tests passing, 0 failing

Run tests:
```bash
pytest                    # All tests
pytest tests/unit/        # Unit tests only
pytest tests/integration/ # Integration tests only
```

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests for new functionality
5. Ensure all tests pass
6. Submit a pull request

## License

This project is for research purposes. Please see LICENSE file for details.

## Research Context

This tool is part of a research project investigating TIM barrel protein structures and their evolutionary relationships. TIM barrels are one of the most common protein folds and are found in many essential enzymes across all domains of life.

## Contact

For questions or collaboration opportunities, please open an issue on GitHub.