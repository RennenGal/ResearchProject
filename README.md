# Protein Data Collector

A comprehensive bioinformatics system for collecting, storing, and analyzing protein data with a focus on TIM barrel proteins. This project provides automated data collection from InterPro and UniProt databases, with support for detailed isoform analysis and structural annotations.

## 🚀 Quick Installation

Get started in minutes with our automated installation scripts:

### Linux/macOS
```bash
git clone https://github.com/RennenGal/ResearchProject.git
cd ResearchProject
./scripts/install.sh
```

### Windows (Command Prompt)
```cmd
git clone https://github.com/RennenGal/ResearchProject.git
cd ResearchProject
scripts\install.bat
```

### Windows (PowerShell)
```powershell
git clone https://github.com/RennenGal/ResearchProject.git
cd ResearchProject
.\scripts\install.ps1
```

The installation scripts will automatically:
- ✅ Check Python 3.8+ and Git installation
- ✅ Create and activate virtual environment
- ✅ Install all dependencies
- ✅ Set up configuration files
- ✅ Run tests to verify installation
- ✅ Provide next steps guidance

For manual installation or troubleshooting, see [docs/prerequisites.md](docs/prerequisites.md).

## ✨ Features

- **🔄 Automated Data Collection**: Collect TIM barrel protein families and human protein data from InterPro and UniProt
- **🗄️ Comprehensive Database Schema**: SQLite database with 67 UniProt fields across 9 categories
- **🧬 Isoform Analysis**: Detailed protein isoform data including sequence, exon annotations, and TIM barrel locations
- **🌐 RESTful API**: Query and export protein data through a FastAPI-based web service
- **🧪 Property-Based Testing**: Comprehensive test suite with 86 tests ensuring data integrity
- **⚡ Rate Limiting**: Built-in rate limiting and retry mechanisms for API compliance
- **⚙️ Flexible Configuration**: Environment-based configuration with support for different deployment scenarios
- **🔗 REST API Only**: Simplified architecture using UniProt REST API for reliable data access

## Quick Start

### Prerequisites

- Python 3.8 or higher
- pip package manager
- SQLite (included with Python)

### Installation

1. **Clone the repository**
   ```bash
   git clone https://github.com/yourusername/protein-data-collector.git
   cd protein-data-collector
   ```

2. **Create a virtual environment**
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. **Install dependencies**
   ```bash
   pip install -r requirements.txt
   ```

4. **Set up configuration**
   ```bash
   cp config.json.example config.json
   cp .env.example .env
   ```
   Edit these files to match your environment (see [Configuration](#configuration) section).

5. **Initialize the database**
   ```bash
   python scripts/create_sqlite_simple.py
   ```

### Basic Usage

1. **Collect TIM barrel entries**
   ```bash
   python scripts/collect_tim_barrel_entries.py
   ```

2. **Collect human proteins**
   ```bash
   python scripts/collect_human_proteins.py
   ```

3. **Start the API server**
   ```bash
   python -m protein_data_collector.server
   ```

4. **Query the data**
   ```bash
   curl http://localhost:8000/api/proteins?limit=10
   ```

## Database Schema

The system uses a three-tier hierarchical structure:

1. **TIM Barrel Entries** (`tim_barrel_entries`): PFAM families and InterPro entries with TIM barrel annotations
2. **InterPro Proteins** (`interpro_proteins`): Human proteins belonging to TIM barrel families
3. **Protein Isoforms** (`proteins`): Detailed isoform data with 67 UniProt fields

### Key Features
- **67 UniProt Fields**: Comprehensive protein data across 9 categories
- **Hierarchical Relationships**: Foreign key constraints maintain data integrity
- **Optimized Indexes**: Efficient querying for common use cases
- **SQLite Backend**: Lightweight, serverless database perfect for research

For detailed schema documentation, see [docs/database-schema.md](docs/database-schema.md).

## API Documentation

The system provides a RESTful API for querying protein data:

### Endpoints

- `GET /api/proteins` - Query protein isoforms with filtering
- `GET /api/proteins/{protein_id}` - Get specific protein details
- `GET /api/families` - List TIM barrel families
- `GET /api/statistics` - Get database statistics
- `GET /health` - Health check endpoint

### Example Queries

```bash
# Get proteins from a specific family
curl "http://localhost:8000/api/proteins?pfam_family=PF00121"

# Get proteins with sequence length filter
curl "http://localhost:8000/api/proteins?min_length=200&max_length=500"

# Get protein statistics
curl "http://localhost:8000/api/statistics"
```

## Configuration

### Environment Variables

Create a `.env` file based on `.env.example`:

```bash
# Database Configuration
DATABASE_TYPE=sqlite
DATABASE_PATH=db/protein_data.db

# API Configuration
INTERPRO_BASE_URL=https://www.ebi.ac.uk/interpro/api/
UNIPROT_BASE_URL=https://rest.uniprot.org/

# Rate Limiting
INTERPRO_REQUESTS_PER_SECOND=10.0
UNIPROT_REQUESTS_PER_SECOND=10.0

# Logging
LOG_LEVEL=INFO
LOG_FORMAT=text
```

### Configuration File

Edit `config.json` for detailed configuration:

```json
{
  "database": {
    "type": "sqlite",
    "path": "db/protein_data.db"
  },
  "api": {
    "interpro_base_url": "https://www.ebi.ac.uk/interpro/api/",
    "uniprot_base_url": "https://rest.uniprot.org/"
  },
  "rate_limiting": {
    "interpro_requests_per_second": 10.0,
    "uniprot_requests_per_second": 10.0
  }
}
```

## Testing

The project includes a comprehensive test suite with both unit tests and property-based tests.

### Run All Tests

```bash
python -m pytest tests/ -v
```

### Run Specific Test Categories

```bash
# API integration tests
python -m pytest tests/test_api_integration.py -v

# Property-based tests
python -m pytest tests/test_complete_isoform_data_collection.py -v

# Database tests
python -m pytest tests/test_database_schema_integrity.py -v
```

### Test Configuration

Tests use a separate configuration file (`config.test.json`) to avoid interfering with production data.

## Data Collection Workflow

1. **TIM Barrel Entries**: Collect PFAM families and InterPro entries with TIM barrel annotations
2. **Human Proteins**: For each TIM barrel entry, collect associated human proteins from InterPro
3. **Protein Isoforms**: For each human protein, collect detailed isoform data from UniProt

### Collection Scripts

- `scripts/collect_tim_barrel_entries.py` - Collect TIM barrel families and InterPro entries
- `scripts/collect_human_proteins.py` - Collect human proteins for TIM barrel families

## Architecture

The system is built with a modular architecture:

- **Data Models** (`protein_data_collector/models/`): Pydantic models for data validation
- **Database Layer** (`protein_data_collector/database/`): SQLAlchemy models and connection management
- **API Clients** (`protein_data_collector/api/`): InterPro and UniProt API clients with rate limiting
- **Collectors** (`protein_data_collector/collector/`): Data collection orchestration
- **Query Engine** (`protein_data_collector/query/`): Database querying and export functionality
- **Web API** (`protein_data_collector/server.py`): FastAPI-based REST API

## Contributing

We welcome contributions! Please see [docs/contributing.md](docs/contributing.md) for guidelines.

### Development Setup

1. Fork the repository
2. Create a feature branch
3. Install development dependencies: `pip install -r requirements-dev.txt`
4. Make your changes
5. Run tests: `python -m pytest`
6. Submit a pull request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- **InterPro Database**: Protein family and domain data
- **UniProt Database**: Comprehensive protein information
- **PFAM Database**: Protein family classifications
- **FastAPI**: Modern web framework for the API
- **SQLAlchemy**: Database ORM and migrations
- **Pydantic**: Data validation and serialization

## Support

For questions, issues, or contributions:

- **Issues**: [GitHub Issues](https://github.com/yourusername/protein-data-collector/issues)
- **Documentation**: [docs/](docs/) directory
- **API Documentation**: Start the server and visit `http://localhost:8000/docs`

## Project Status

This project is actively maintained and used for TIM barrel protein research. The database currently contains:

- 49 TIM barrel entries (18 PFAM families + 31 InterPro entries)
- 407 human proteins associated with TIM barrel families
- Comprehensive isoform data with structural annotations

For the latest statistics, run the collection scripts or query the `/api/statistics` endpoint.