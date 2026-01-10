# Project Organization Guidelines

## Directory Structure Standards

### Documentation (`docs/`)
- **Purpose**: Long-term, permanent documentation only
- **Contents**: 
  - API interfaces and specifications
  - Database schema documentation
  - User guides and workflows
  - Architecture documentation
- **Exclusions**: Temporary analysis, migration summaries, implementation guides

### Temporary Documentation (`docs/temp/`)
- **Purpose**: Temporary documentation, analysis, and implementation notes
- **Contents**:
  - Migration summaries and logs
  - Implementation guides and analysis
  - Project completion summaries
  - Temporary research and exploration docs
- **Lifecycle**: Review periodically and either promote to `docs/` or archive/delete

### Scripts (`scripts/`)
- **Purpose**: Production-ready, reusable scripts only
- **Contents**:
  - Data collection scripts
  - Database management utilities
  - Core functionality scripts
- **Quality**: Well-documented, error-handled, production-ready

### Temporary Scripts (`scripts/temp/`)
- **Purpose**: Development, testing, and one-time use scripts
- **Contents**:
  - Debug and testing scripts
  - One-time migration scripts
  - Experimental code
  - API testing and exploration
- **Lifecycle**: Review periodically and either promote to `scripts/` or delete

## File Organization Rules

### Documentation Files
1. **Permanent docs** → `docs/`
   - `database-schema.md`
   - `interpro-api-interface.md`
   - `scripts-workflow-guide.md`

2. **Temporary docs** → `docs/temp/`
   - Migration summaries
   - Implementation guides
   - Analysis documents
   - Project completion summaries

### Script Files
1. **Production scripts** → `scripts/`
   - `collect_tim_barrel_entries.py`
   - `collect_human_proteins.py`
   - Database creation and management scripts

2. **Temporary scripts** → `scripts/temp/`
   - Debug and test scripts
   - One-time migration scripts
   - API exploration scripts

## Maintenance Guidelines

### Regular Cleanup
- **Monthly**: Review `docs/temp/` and `scripts/temp/`
- **Decision criteria**:
  - Still needed? Keep in temp
  - Permanent value? Move to main directory
  - Obsolete? Delete

### File Naming
- Use descriptive, consistent naming
- Include purpose in temporary file names
- Date temporary files when relevant

### Documentation Standards
- All permanent docs must be well-structured
- Include purpose and scope at the top
- Maintain table of contents for longer docs
- Use consistent markdown formatting

## Implementation Notes
- This structure supports clean project maintenance
- Separates permanent assets from temporary work
- Enables easy cleanup and organization
- Follows software engineering best practices