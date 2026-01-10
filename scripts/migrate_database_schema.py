#!/usr/bin/env python3
"""
Database migration script to update schema from old structure to new unified structure.

This script:
1. Backs up existing data from pfam_families and interpro_proteins tables
2. Creates new unified tim_barrel_entries table
3. Migrates data to new structure with composite primary keys
4. Updates interpro_proteins table to use composite primary key
5. Handles foreign key constraints properly

Usage:
    python scripts/migrate_database_schema.py [--config CONFIG_FILE] [--dry-run]
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import List, Dict, Any
from datetime import datetime

# Add the project root to Python path
sys.path.insert(0, str(Path(__file__).parent.parent))

from protein_data_collector.config import load_config_from_file, set_config, LoggingConfig
from protein_data_collector.database.connection import get_database_manager
from protein_data_collector.logging_config import setup_logging
from sqlalchemy import text


class DatabaseMigrator:
    """Handles database schema migration from old to new structure."""
    
    def __init__(self, dry_run: bool = False):
        """
        Initialize the migrator.
        
        Args:
            dry_run: If True, don't actually modify database
        """
        self.dry_run = dry_run
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
        self.db_manager = get_database_manager()
        
        # Migration statistics
        self.stats = {
            'pfam_families_migrated': 0,
            'interpro_entries_created': 0,
            'proteins_migrated': 0,
            'errors': 0
        }
    
    def check_current_schema(self) -> Dict[str, bool]:
        """
        Check what tables currently exist in the database.
        
        Returns:
            Dictionary indicating which tables exist
        """
        self.logger.info("Checking current database schema")
        
        with self.db_manager.get_session() as session:
            # Check for existing tables
            result = session.execute(text("""
                SELECT table_name 
                FROM information_schema.tables 
                WHERE table_schema = DATABASE()
            """))
            
            existing_tables = {row[0] for row in result.fetchall()}
            
            schema_status = {
                'has_pfam_families': 'pfam_families' in existing_tables,
                'has_tim_barrel_entries': 'tim_barrel_entries' in existing_tables,
                'has_old_interpro_proteins': False,
                'has_new_interpro_proteins': False
            }
            
            # Check interpro_proteins table structure
            if 'interpro_proteins' in existing_tables:
                # Check if it has the old structure (single primary key)
                result = session.execute(text("""
                    SELECT COLUMN_NAME, COLUMN_KEY
                    FROM information_schema.columns 
                    WHERE table_schema = DATABASE() 
                    AND table_name = 'interpro_proteins'
                    AND COLUMN_KEY = 'PRI'
                """))
                
                primary_keys = [row[0] for row in result.fetchall()]
                
                if len(primary_keys) == 1 and primary_keys[0] == 'uniprot_id':
                    schema_status['has_old_interpro_proteins'] = True
                elif len(primary_keys) == 2 and 'uniprot_id' in primary_keys and 'tim_barrel_accession' in primary_keys:
                    schema_status['has_new_interpro_proteins'] = True
            
            self.logger.info(f"Schema status: {schema_status}")
            return schema_status
    
    def backup_existing_data(self) -> Dict[str, List[Dict[str, Any]]]:
        """
        Backup existing data before migration.
        
        Returns:
            Dictionary containing backed up data
        """
        self.logger.info("Backing up existing data")
        
        backup_data = {
            'pfam_families': [],
            'interpro_proteins': []
        }
        
        with self.db_manager.get_session() as session:
            # Backup PFAM families
            try:
                result = session.execute(text("SELECT * FROM pfam_families"))
                columns = result.keys()
                backup_data['pfam_families'] = [
                    dict(zip(columns, row)) for row in result.fetchall()
                ]
                self.logger.info(f"Backed up {len(backup_data['pfam_families'])} PFAM families")
            except Exception as e:
                self.logger.warning(f"Could not backup pfam_families: {e}")
            
            # Backup InterPro proteins
            try:
                result = session.execute(text("SELECT * FROM interpro_proteins"))
                columns = result.keys()
                backup_data['interpro_proteins'] = [
                    dict(zip(columns, row)) for row in result.fetchall()
                ]
                self.logger.info(f"Backed up {len(backup_data['interpro_proteins'])} InterPro proteins")
            except Exception as e:
                self.logger.warning(f"Could not backup interpro_proteins: {e}")
        
        return backup_data
    
    def create_new_schema(self) -> None:
        """Create the new unified schema structure."""
        self.logger.info("Creating new unified schema structure")
        
        if self.dry_run:
            self.logger.info("DRY RUN: Would create new schema structure")
            return
        
        with self.db_manager.get_session() as session:
            # Create tim_barrel_entries table
            session.execute(text("""
                CREATE TABLE IF NOT EXISTS tim_barrel_entries (
                    accession VARCHAR(20) PRIMARY KEY,
                    entry_type VARCHAR(20) NOT NULL,
                    name VARCHAR(255) NOT NULL,
                    description TEXT,
                    interpro_type VARCHAR(50),
                    tim_barrel_annotation TEXT NOT NULL,
                    member_databases JSON,
                    interpro_id VARCHAR(20),
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    
                    INDEX idx_tim_barrel_type (entry_type),
                    INDEX idx_tim_barrel_name (name),
                    INDEX idx_tim_barrel_interpro (interpro_id),
                    INDEX idx_tim_barrel_created (created_at)
                ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci
            """))
            
            session.commit()
            self.logger.info("Created tim_barrel_entries table")
    
    def migrate_pfam_families(self, backup_data: Dict[str, List[Dict[str, Any]]]) -> None:
        """
        Migrate PFAM families to the new tim_barrel_entries table.
        
        Args:
            backup_data: Backed up data from old tables
        """
        self.logger.info("Migrating PFAM families to tim_barrel_entries")
        
        pfam_families = backup_data.get('pfam_families', [])
        
        if not pfam_families:
            self.logger.info("No PFAM families to migrate")
            return
        
        if self.dry_run:
            self.logger.info(f"DRY RUN: Would migrate {len(pfam_families)} PFAM families")
            self.stats['pfam_families_migrated'] = len(pfam_families)
            return
        
        with self.db_manager.get_session() as session:
            for family in pfam_families:
                try:
                    # Insert into tim_barrel_entries with entry_type = 'pfam'
                    session.execute(text("""
                        INSERT INTO tim_barrel_entries 
                        (accession, entry_type, name, description, tim_barrel_annotation, interpro_id, created_at)
                        VALUES (:accession, :entry_type, :name, :description, :tim_barrel_annotation, :interpro_id, :created_at)
                        ON DUPLICATE KEY UPDATE
                        name = VALUES(name),
                        description = VALUES(description),
                        tim_barrel_annotation = VALUES(tim_barrel_annotation),
                        interpro_id = VALUES(interpro_id)
                    """), {
                        'accession': family['accession'],
                        'entry_type': 'pfam',
                        'name': family['name'],
                        'description': family.get('description'),
                        'tim_barrel_annotation': family['tim_barrel_annotation'],
                        'interpro_id': family.get('interpro_id'),
                        'created_at': family.get('created_at', datetime.now())
                    })
                    
                    self.stats['pfam_families_migrated'] += 1
                    
                except Exception as e:
                    self.logger.error(f"Failed to migrate PFAM family {family.get('accession')}: {e}")
                    self.stats['errors'] += 1
            
            session.commit()
            self.logger.info(f"Migrated {self.stats['pfam_families_migrated']} PFAM families")
    
    def recreate_interpro_proteins_table(self, backup_data: Dict[str, List[Dict[str, Any]]]) -> None:
        """
        Recreate interpro_proteins table with composite primary key.
        
        Args:
            backup_data: Backed up data from old tables
        """
        self.logger.info("Recreating interpro_proteins table with composite primary key")
        
        proteins = backup_data.get('interpro_proteins', [])
        
        if self.dry_run:
            self.logger.info(f"DRY RUN: Would recreate interpro_proteins table and migrate {len(proteins)} proteins")
            self.stats['proteins_migrated'] = len(proteins)
            return
        
        with self.db_manager.get_session() as session:
            # Drop existing interpro_proteins table (this will also drop dependent tables)
            self.logger.warning("Dropping existing interpro_proteins table and dependent data")
            session.execute(text("SET FOREIGN_KEY_CHECKS = 0"))
            session.execute(text("DROP TABLE IF EXISTS proteins"))
            session.execute(text("DROP TABLE IF EXISTS interpro_proteins"))
            session.execute(text("SET FOREIGN_KEY_CHECKS = 1"))
            
            # Create new interpro_proteins table with composite primary key
            session.execute(text("""
                CREATE TABLE interpro_proteins (
                    uniprot_id VARCHAR(20) NOT NULL,
                    tim_barrel_accession VARCHAR(20) NOT NULL,
                    name VARCHAR(255),
                    organism VARCHAR(100) DEFAULT 'Homo sapiens',
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    
                    PRIMARY KEY (uniprot_id, tim_barrel_accession),
                    FOREIGN KEY (tim_barrel_accession) REFERENCES tim_barrel_entries(accession) ON DELETE CASCADE,
                    INDEX idx_interpro_tim_barrel (tim_barrel_accession),
                    INDEX idx_interpro_organism (organism),
                    INDEX idx_interpro_name (name),
                    INDEX idx_interpro_created (created_at)
                ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci
            """))
            
            # Migrate protein data
            for protein in proteins:
                try:
                    # Map old pfam_accession to tim_barrel_accession
                    tim_barrel_accession = protein.get('pfam_accession')
                    if not tim_barrel_accession:
                        self.logger.warning(f"Protein {protein.get('uniprot_id')} has no pfam_accession, skipping")
                        continue
                    
                    session.execute(text("""
                        INSERT INTO interpro_proteins 
                        (uniprot_id, tim_barrel_accession, name, organism, created_at)
                        VALUES (:uniprot_id, :tim_barrel_accession, :name, :organism, :created_at)
                        ON DUPLICATE KEY UPDATE
                        name = VALUES(name),
                        organism = VALUES(organism)
                    """), {
                        'uniprot_id': protein['uniprot_id'],
                        'tim_barrel_accession': tim_barrel_accession,
                        'name': protein.get('name'),
                        'organism': protein.get('organism', 'Homo sapiens'),
                        'created_at': protein.get('created_at', datetime.now())
                    })
                    
                    self.stats['proteins_migrated'] += 1
                    
                except Exception as e:
                    self.logger.error(f"Failed to migrate protein {protein.get('uniprot_id')}: {e}")
                    self.stats['errors'] += 1
            
            session.commit()
            self.logger.info(f"Migrated {self.stats['proteins_migrated']} proteins to new structure")
    
    def cleanup_old_tables(self) -> None:
        """Remove old tables after successful migration."""
        self.logger.info("Cleaning up old tables")
        
        if self.dry_run:
            self.logger.info("DRY RUN: Would drop old pfam_families table")
            return
        
        with self.db_manager.get_session() as session:
            # Drop old pfam_families table
            session.execute(text("DROP TABLE IF EXISTS pfam_families"))
            session.commit()
            self.logger.info("Dropped old pfam_families table")
    
    def run_migration(self) -> None:
        """Run the complete database migration."""
        self.logger.info("Starting database schema migration")
        
        try:
            # Check current schema
            schema_status = self.check_current_schema()
            
            # If already migrated, skip
            if schema_status['has_tim_barrel_entries'] and schema_status['has_new_interpro_proteins']:
                self.logger.info("Database already migrated to new schema")
                return
            
            # If no old data, just create new schema
            if not schema_status['has_pfam_families'] and not schema_status['has_old_interpro_proteins']:
                self.logger.info("No existing data found, creating fresh schema")
                self.create_new_schema()
                return
            
            # Backup existing data
            backup_data = self.backup_existing_data()
            
            # Create new schema
            self.create_new_schema()
            
            # Migrate PFAM families
            self.migrate_pfam_families(backup_data)
            
            # Recreate interpro_proteins table
            self.recreate_interpro_proteins_table(backup_data)
            
            # Cleanup old tables
            self.cleanup_old_tables()
            
            self.logger.info("Database migration completed successfully")
            
        except Exception as e:
            self.logger.error(f"Migration failed: {e}")
            raise
    
    def print_summary(self) -> None:
        """Print migration summary statistics."""
        print("\n" + "="*70)
        print("DATABASE MIGRATION SUMMARY")
        print("="*70)
        print(f"PFAM families migrated:            {self.stats['pfam_families_migrated']}")
        print(f"Proteins migrated:                 {self.stats['proteins_migrated']}")
        print(f"Total errors:                      {self.stats['errors']}")
        print("="*70)
        
        if self.stats['errors'] > 0:
            print(f"⚠️  {self.stats['errors']} errors occurred during migration")
        else:
            print("✅ Migration completed successfully!")
        
        if self.dry_run:
            print("🔍 DRY RUN MODE - No data was actually modified")


def main():
    """Main script entry point."""
    parser = argparse.ArgumentParser(
        description="Migrate database schema from old to new unified structure"
    )
    parser.add_argument(
        '--config', 
        type=str, 
        default='config.test.json',
        help='Configuration file path (default: config.test.json)'
    )
    parser.add_argument(
        '--dry-run', 
        action='store_true',
        help='Run without modifying database'
    )
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose logging'
    )
    
    args = parser.parse_args()
    
    # Setup logging
    log_level = 'DEBUG' if args.verbose else 'INFO'
    log_config = LoggingConfig(
        level=log_level,
        format='text',
        log_file=None
    )
    setup_logging(log_config)
    
    logger = logging.getLogger(__name__)
    
    try:
        # Load configuration
        config_path = Path(args.config)
        if not config_path.exists():
            logger.error(f"Configuration file not found: {config_path}")
            sys.exit(1)
        
        logger.info(f"Loading configuration from: {config_path}")
        config = load_config_from_file(str(config_path))
        set_config(config)
        
        # Test database connection
        db_manager = get_database_manager()
        if not db_manager.test_connection():
            logger.error("Database connection test failed")
            sys.exit(1)
        
        logger.info("Database connection successful")
        
        # Create migrator and run migration
        migrator = DatabaseMigrator(dry_run=args.dry_run)
        migrator.run_migration()
        
        # Print summary
        migrator.print_summary()
        
    except KeyboardInterrupt:
        logger.info("Migration interrupted by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Migration failed: {str(e)}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()