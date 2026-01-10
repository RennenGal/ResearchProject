# Design Document: Protein Data Collector

## Overview

The Protein Data Collector is a bioinformatics application that systematically collects and analyzes TIM barrel proteins from humans using a hierarchical data collection approach. The system integrates with InterPro and UniProt databases to build a comprehensive local database of protein families, proteins, and their isoforms, enabling efficient querying and analysis of TIM barrel structural annotations.

The system follows a three-tier data collection hierarchy:
1. **PFAM Families** - Identify families with TIM barrel annotations from InterPro
2. **Human Proteins** - Find proteins belonging to those families (inheriting TIM barrel annotations)
3. **Protein Isoforms** - Collect detailed sequence and structural data from UniProt

## Architecture

The system employs a modular, service-oriented architecture with clear separation of concerns:

```mermaid
graph TB
    subgraph "External APIs"
        InterPro[InterPro API]
        UniProt[UniProt REST API]
    end
    
    subgraph "Core Services"
        DC[Data Collector]
        VA[Validation Service]
        RC[Retry Controller]
    end
    
    subgraph "Data Layer"
        DB[(Local Database)]
        Cache[Response Cache]
    end
    
    subgraph "Query Interface"
        QE[Query Engine]
        API[REST API]
        CLI[CLI Interface]
    end
    
    InterPro --> DC
    UniProt --> DC
    DC --> VA
    DC --> RC
    VA --> DB
    DB --> QE
    QE --> API
    QE --> CLI
    DC --> Cache
    Cache --> DC
```

### Key Architectural Principles

- **Hierarchical Data Flow**: Data collection follows the PFAM → Protein → Isoform hierarchy
- **Fault Tolerance**: Configurable retry mechanisms with exponential backoff for all external API calls
- **Data Integrity**: Comprehensive validation at each collection stage
- **Performance**: Local caching and efficient database design for fast queries
- **Extensibility**: Modular design supports adding new data sources and analysis capabilities

## Components and Interfaces

### 1. Data Collector Service

**Responsibilities:**
- Orchestrate the three-phase data collection process
- Manage API interactions with InterPro and UniProt REST API
- Handle rate limiting and retry logic
- Coordinate data validation and storage

**Key Interfaces:**
```python
class DataCollector:
    def collect_pfam_families(self, tim_barrel_query: str) -> List[PfamFamily]
    def collect_human_proteins(self, pfam_families: List[PfamFamily]) -> List[InterProProtein]
    def collect_protein_isoforms(self, interpro_proteins: List[InterProProtein]) -> List[Protein]
    def run_full_collection(self) -> CollectionReport
```

### 2. API Integration Layer

**InterPro Integration:**
- Base URL: `https://www.ebi.ac.uk/interpro/api/`
- Primary endpoints:
  - `/entry/pfam/` - List PFAM families
  - `/entry/pfam/{accession}/` - Get family details
  - `/protein/UniProt/entry/pfam/{accession}/` - Get proteins in family
- Query parameters for TIM barrel annotation filtering

**UniProt Integration:**
- Direct REST API approach: `https://rest.uniprot.org/uniprotkb/`
- Key operations: protein lookup, isoform retrieval, sequence/annotation extraction

**Unified API Interface:**
```python
class APIClient:
    def query_interpro(self, endpoint: str, params: dict) -> dict
    def query_uniprot_rest(self, accession: str, fields: List[str]) -> dict
    def with_retry(self, operation: Callable, max_retries: int) -> Any
```

### 3. Retry Controller

**Configuration:**
- Configurable retry count K (default: 3)
- Exponential backoff: initial_delay * (backoff_multiplier ^ attempt)
- Maximum delay cap to prevent excessive wait times
- Per-database retry policies

**Implementation:**
```python
@dataclass
class RetryConfig:
    max_retries: int = 3
    initial_delay: float = 1.0
    backoff_multiplier: float = 2.0
    max_delay: float = 60.0
    
class RetryController:
    def execute_with_retry(self, operation: Callable, config: RetryConfig) -> Any
    def log_retry_attempt(self, database: str, operation: str, attempt: int)
```

### 4. Validation Service

**Validation Rules:**
- Protein sequences: Valid amino acid characters only (ACDEFGHIKLMNPQRSTVWY)
- TIM barrel coordinates: Within protein sequence bounds
- Data completeness: Required fields present and non-empty
- Cross-reference integrity: Valid relationships between entities

**Interface:**
```python
class ValidationService:
    def validate_protein_sequence(self, sequence: str) -> ValidationResult
    def validate_tim_barrel_location(self, location: dict, sequence_length: int) -> ValidationResult
    def validate_data_completeness(self, entity: dict, required_fields: List[str]) -> ValidationResult
```

### 5. Local Database

**Schema Design:**
```sql
-- PFAM families with TIM barrel annotations
CREATE TABLE pfam_families (
    accession VARCHAR(20) PRIMARY KEY,
    name VARCHAR(255) NOT NULL,
    description TEXT,
    tim_barrel_annotation TEXT NOT NULL,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- Human proteins belonging to PFAM families (from InterPro)
CREATE TABLE interpro_proteins (
    uniprot_id VARCHAR(20) PRIMARY KEY,
    pfam_accession VARCHAR(20) NOT NULL,
    name VARCHAR(255),
    organism VARCHAR(100) DEFAULT 'Homo sapiens',
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (pfam_accession) REFERENCES pfam_families(accession)
);

-- Protein isoforms with detailed annotations (from UniProt)
CREATE TABLE proteins (
    isoform_id VARCHAR(30) PRIMARY KEY,
    parent_protein_id VARCHAR(20) NOT NULL,
    sequence TEXT NOT NULL,
    sequence_length INTEGER NOT NULL,
    exon_annotations JSON,
    exon_count INTEGER,
    tim_barrel_location JSON,
    organism VARCHAR(100),
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (parent_protein_id) REFERENCES interpro_proteins(uniprot_id)
);

-- Indexes for efficient querying
CREATE INDEX idx_interpro_proteins_pfam ON interpro_proteins(pfam_accession);
CREATE INDEX idx_proteins_parent ON proteins(parent_protein_id);
CREATE INDEX idx_interpro_proteins_organism ON interpro_proteins(organism);
```

### 6. Query Engine

**Query Capabilities:**
- Search by PFAM family identifier
- Filter by TIM barrel structural features
- Retrieve protein isoforms
- Complex filtering (sequence length, exon count, organism)
- Data quality metrics

**Interface:**
```python
class QueryEngine:
    def search_by_pfam_family(self, pfam_id: str) -> List[Protein]
    def search_by_tim_barrel_features(self, criteria: dict) -> List[Protein]
    def get_protein_isoforms(self, protein_id: str) -> List[ProteinIsoform]
    def filter_proteins(self, filters: dict) -> List[Protein]
    def export_results(self, results: List[Any], format: str) -> str
```

## Data Models

### Core Entities

```python
@dataclass
class PfamFamily:
    accession: str
    name: str
    description: str
    tim_barrel_annotation: str
    interpro_id: Optional[str] = None
    created_at: datetime = field(default_factory=datetime.now)

@dataclass
class InterProProtein:
    uniprot_id: str
    pfam_accession: str
    name: str
    organism: str = "Homo sapiens"
    basic_metadata: dict = field(default_factory=dict)
    created_at: datetime = field(default_factory=datetime.now)

@dataclass
class Protein:
    isoform_id: str
    parent_protein_id: str
    sequence: str
    sequence_length: int
    exon_annotations: dict
    exon_count: int
    tim_barrel_location: dict
    organism: str
    name: str
    description: str
    created_at: datetime = field(default_factory=datetime.now)

@dataclass
class CollectionReport:
    pfam_families_found: int
    interpro_proteins_collected: int
    protein_isoforms_collected: int
    validation_errors: List[str]
    collection_duration: timedelta
    success_rate: float
```

## Error Handling

### Error Categories

1. **Network Errors**: Connection timeouts, DNS failures
2. **API Errors**: Rate limiting (429), server errors (5xx), authentication issues
3. **Data Errors**: Invalid responses, missing required fields, malformed data
4. **Validation Errors**: Invalid sequences, coordinate mismatches, constraint violations

### Error Handling Strategy

```python
class ErrorHandler:
    def handle_network_error(self, error: NetworkError) -> ErrorAction
    def handle_api_error(self, error: APIError) -> ErrorAction
    def handle_data_error(self, error: DataError) -> ErrorAction
    def handle_validation_error(self, error: ValidationError) -> ErrorAction

@enum
class ErrorAction:
    RETRY = "retry"
    SKIP = "skip"
    FAIL = "fail"
    LOG_AND_CONTINUE = "log_and_continue"
```

### Logging Strategy

- **Structured logging** with JSON format for machine readability
- **Log levels**: DEBUG (API calls), INFO (progress), WARN (recoverable errors), ERROR (failures)
- **Contextual information**: operation type, database, entity IDs, retry attempts
- **Performance metrics**: response times, success rates, collection statistics

## Testing Strategy

The testing approach combines unit tests for specific functionality with property-based tests for comprehensive validation of system behavior across diverse inputs.

### Unit Testing

**Focus Areas:**
- API client functionality with mocked responses
- Data validation logic with known good/bad inputs
- Database operations with test fixtures
- Error handling with simulated failure conditions
- Configuration management and retry logic

**Key Test Cases:**
- InterPro API response parsing
- UniProt REST API integration
- Protein sequence validation
- TIM barrel coordinate validation
- Database schema constraints
- Retry mechanism behavior

### Property-Based Testing

Property-based tests will validate universal properties that should hold across all valid system inputs, using a minimum of 100 iterations per test to ensure comprehensive coverage through randomization.

**Testing Framework**: Hypothesis (Python) - provides powerful generators for biological data and supports complex property definitions.

**Test Configuration**: Each property test will be tagged with format: `Feature: protein-data-collector, Property {number}: {property_text}` and reference the corresponding design document property.

### Integration Testing

**End-to-End Scenarios:**
- Complete data collection workflow from PFAM discovery to isoform storage
- Database integrity during concurrent operations
- Query engine performance with realistic data volumes

**Test Data Management:**
- Curated test datasets with known TIM barrel proteins
- Mock API responses for consistent testing
- Database fixtures for query testing
- Performance benchmarks for regression testing

## Correctness Properties

*A property is a characteristic or behavior that should hold true across all valid executions of a system—essentially, a formal statement about what the system should do. Properties serve as the bridge between human-readable specifications and machine-verifiable correctness guarantees.*

The following properties define the correctness criteria for the Protein Data Collector system. Each property represents a universal rule that must hold across all valid inputs and system states.

### Property 1: PFAM Family Collection Completeness
*For any* InterPro query response containing PFAM families with TIM barrel annotations, the Data_Collector should retrieve exactly those families that contain TIM barrel structural annotations and no others.
**Validates: Requirements 1.1, 1.3**

### Property 2: Hierarchical Data Integrity
*For any* collected dataset, every protein should be correctly linked to exactly one PFAM family, and every isoform should be correctly linked to exactly one protein, maintaining the three-tier hierarchy.
**Validates: Requirements 2.3, 3.7, 4.4**

### Property 3: Human Organism Filtering
*For any* protein collection operation, all retrieved proteins should have organism field set to "Homo sapiens" and no proteins from other organisms should be included.
**Validates: Requirements 2.2**

### Property 4: Complete Isoform Data Collection
*For any* protein isoform retrieved from UniProt, the system should collect all required fields: sequence information, exon annotations, TIM barrel location, protein name, description, and organism information.
**Validates: Requirements 3.2, 3.3, 3.4, 3.5, 3.6**

### Property 5: Configurable Retry Behavior
*For any* external API call (InterPro, UniProt REST API) that fails, the system should retry exactly K times (where K is configurable) using exponential backoff, then log the failure and continue processing.
**Validates: Requirements 1.4, 3.8, 9.1, 9.2, 9.5**

### Property 6: Data Validation Consistency
*For any* protein sequence, the system should accept only sequences containing valid amino acid characters (ACDEFGHIKLMNPQRSTVWY), and for any TIM barrel location coordinates, they should be within the bounds of the protein sequence length.
**Validates: Requirements 7.1, 7.2**

### Property 7: Database Storage Completeness
*For any* successfully validated entity (PFAM family, protein, or isoform), the system should store all required fields in the local database with proper relational integrity maintained.
**Validates: Requirements 4.1, 4.2, 4.3, 4.5**

### Property 8: Query Result Accuracy
*For any* query by PFAM family identifier, TIM barrel features, or protein identifier, the Query_Engine should return exactly the entities that match the specified criteria and include all required display fields.
**Validates: Requirements 5.1, 5.2, 5.3, 5.5**

### Property 9: UniProt REST API Reliability
*For any* UniProt data request, the system should use the REST API with comprehensive error handling and retry logic to ensure reliable data retrieval.
**Validates: Requirements 8.1, 8.2, 8.3**

### Property 10: Rate Limiting Compliance
*For any* sequence of API calls to external databases, the system should respect rate limiting policies by spacing requests appropriately and implementing exponential backoff when rate limits are exceeded.
**Validates: Requirements 6.1, 6.2**

### Property 11: Data Export Format Validity
*For any* query result set, when exported in standard bioinformatics formats (FASTA, JSON, CSV), the output should be valid according to the format specification and contain all requested data fields.
**Validates: Requirements 5.6**

### Property 12: Collection Progress Persistence
*For any* interrupted data collection operation, the system should be able to resume from the last successfully completed step without duplicating work or losing progress.
**Validates: Requirements 6.5**
## Implementation Notes

### Technology Stack Recommendations

**Backend Framework**: Python with FastAPI for REST API endpoints
- Excellent bioinformatics ecosystem (BioPython, requests, pandas)
- Strong async support for concurrent API calls
- Built-in data validation with Pydantic models

**Database**: MySQL with JSON support
- ACID compliance for data integrity
- JSON columns for flexible annotation storage (MySQL 5.7+)
- Excellent performance for complex queries
- Wide adoption and strong ecosystem support

**Testing Framework**: Pytest with Hypothesis for property-based testing
- Mature property-based testing library
- Excellent biological data generators
- Integration with standard Python testing ecosystem

**API Client Libraries**:
- `httpx` for async HTTP requests
- `tenacity` for retry logic implementation
- `pydantic` for data validation and serialization

### Configuration Management

```python
@dataclass
class SystemConfig:
    # Database configuration
    database_url: str  # e.g., "mysql://user:password@localhost/protein_db"
    
    # API configuration
    interpro_base_url: str = "https://www.ebi.ac.uk/interpro/api/"
    uniprot_base_url: str = "https://rest.uniprot.org/"
    
    # Retry configuration
    max_retries: int = 3
    initial_delay: float = 1.0
    backoff_multiplier: float = 2.0
    max_delay: float = 60.0
    
    # Rate limiting
    interpro_requests_per_second: float = 10.0
    uniprot_requests_per_second: float = 5.0
    
    # Collection configuration
    batch_size: int = 100
    cache_ttl_hours: int = 24
```

### Performance Considerations

**Concurrent Processing**: Use asyncio for concurrent API calls while respecting rate limits
**Caching Strategy**: Implement response caching with TTL to reduce API load
**Database Optimization**: Use connection pooling and prepared statements
**Memory Management**: Process data in batches to handle large datasets
**Progress Tracking**: Implement checkpointing for long-running operations

### Security Considerations

**API Key Management**: Store API keys in environment variables or secure key management
**Input Validation**: Validate all external API responses before processing
**SQL Injection Prevention**: Use parameterized queries and ORM
**Rate Limit Compliance**: Respect external API terms of service
**Data Privacy**: Handle protein data according to relevant data protection regulations

### Monitoring and Observability

**Metrics Collection**:
- API response times and success rates
- Database query performance
- Collection progress and throughput
- Error rates by category

**Health Checks**:
- Database connectivity
- External API availability
- Disk space and memory usage

**Alerting**:
- Collection failures
- API rate limit violations
- Database constraint violations
- Unusual error patterns