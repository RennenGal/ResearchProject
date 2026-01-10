#!/usr/bin/env python3
"""
Database cleanup script to remove redundant tables and ensure correct schema.

This script:
1. Checks current database state
2. Removes redundant pfam_families and interpro_entries tables
3. Creates the correct interpro_proteins table with composite primary key
4. Creates the proteins table with proper foreign key relationships

Usage:
    python scripts/cleanup_database_schema.py [--config CONFIG_FILE] [--dry-run]
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


class DatabaseCleaner:
    """Cleans up redundant database tables and ensures correct schema."""
    
    def __init__(self, dry_run: bool = False):
        """
        Initialize the cleaner.
        
        Args:
            dry_run: If True, don't actually modify database
        """
        self.dry_run = dry_run
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
        self.db_manager = get_database_manager()
        
        # Cleanup statistics
        self.stats = {
            'tables_dropped': 0,
            'tables_created': 0,
            'errors': 0
        }
    
    def check_current_tables(self) -> List[str]:
        """
        Check what tables currently exist in the database.
        
        Returns:
            List of existing table names
        """
        self.logger.info("Checking current database tables")
        
        with self.db_manager.get_session() as session:
            result = session.execute(text("SHOW TABLES"))
            tables = [row[0] for row in result.fetchall()]
            
            self.logger.info(f"Found tables: {tables}")
            return tables
    
    def drop_redundant_tables(self, existing_tables: List[str]) -> None:
        """
        Drop redundant tables that are replaced by tim_barrel_entries.
        
        Args:
            existing_tables: List of existing table names
        """
        redundant_tables = ['pfam_families', 'interpro_entries']
        
        if self.dry_run:
            tables_to_drop = [t for t in redundant_tables if t in existing_tables]
            self.logger.info(f"DRY RUN: Would drop redundant tables: {tables_to_drop}")
            self.stats['tables_dropped'] = len(tables_to_drop)
            return
        
        with self.db_manager.get_session() as session:
            # Disable foreign key checks temporarily
            session.execute(text("SET FOREIGN_KEY_CHECKS = 0"))
            
            for table in redundant_tables:
                if table in existing_tables:
                    try:
                        self.logger.info(f"Dropping redundant table: {table}")
                        session.execute(text(f"DROP TABLE IF EXISTS {table}"))
                        self.stats['tables_dropped'] += 1
                    except Exception as e:
                        self.logger.error(f"Failed to drop table {table}: {e}")
                        self.stats['errors'] += 1
            
            # Re-enable foreign key checks
            session.execute(text("SET FOREIGN_KEY_CHECKS = 1"))
            session.commit()
            
            self.logger.info(f"Dropped {self.stats['tables_dropped']} redundant tables")
    
    def create_interpro_proteins_table(self) -> None:
        """Create the interpro_proteins table with composite primary key."""
        self.logger.info("Creating interpro_proteins table with composite primary key")
        
        if self.dry_run:
            self.logger.info("DRY RUN: Would create interpro_proteins table")
            self.stats['tables_created'] += 1
            return
        
        with self.db_manager.get_session() as session:
            try:
                # Create interpro_proteins table with composite primary key
                session.execute(text("""
                    CREATE TABLE IF NOT EXISTS interpro_proteins (
                        uniprot_id VARCHAR(20) NOT NULL,
                        tim_barrel_accession VARCHAR(20) NOT NULL,
                        name VARCHAR(255),
                        organism VARCHAR(100) DEFAULT 'Homo sapiens',
                        created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
                        
                        PRIMARY KEY (uniprot_id, tim_barrel_accession),
                        FOREIGN KEY (tim_barrel_accession) REFERENCES tim_barrel_entries(accession) ON DELETE CASCADE,
                        INDEX idx_interpro_tim_barrel (tim_barrel_accession),
                        INDEX idx_interpro_organism (organism),
                        INDEX idx_interpro_name (name),
                        INDEX idx_interpro_created (created_at)
                    ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci
                """))
                
                session.commit()
                self.stats['tables_created'] += 1
                self.logger.info("Created interpro_proteins table")
                
            except Exception as e:
                self.logger.error(f"Failed to create interpro_proteins table: {e}")
                self.stats['errors'] += 1
                raise
    
    def create_proteins_table(self) -> None:
        """Create the proteins table with proper foreign key relationships."""
        self.logger.info("Creating proteins table with composite foreign key")
        
        if self.dry_run:
            self.logger.info("DRY RUN: Would create proteins table")
            self.stats['tables_created'] += 1
            return
        
        with self.db_manager.get_session() as session:
            try:
                # Create proteins table with composite foreign key
                session.execute(text("""
                    CREATE TABLE IF NOT EXISTS proteins (
                        isoform_id VARCHAR(30) PRIMARY KEY,
                        parent_protein_id VARCHAR(20) NOT NULL,
                        parent_tim_barrel_accession VARCHAR(20) NOT NULL,
                        sequence TEXT NOT NULL,
                        sequence_length INTEGER NOT NULL,
                        exon_annotations JSON,
                        exon_count INTEGER,
                        tim_barrel_location JSON,
                        organism VARCHAR(100),
                        name VARCHAR(255),
                        description TEXT,
                        created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
                        
                        FOREIGN KEY (parent_protein_id, parent_tim_barrel_accession) 
                            REFERENCES interpro_proteins(uniprot_id, tim_barrel_accession) ON DELETE CASCADE,
                        INDEX idx_proteins_parent (parent_protein_id),
                        INDEX idx_proteins_parent_composite (parent_protein_id, parent_tim_barrel_accession),
                        INDEX idx_proteins_organism (organism),
                        INDEX idx_proteins_length (sequence_length),
                        INDEX idx_proteins_exon_count (exon_count),
                        INDEX idx_proteins_name (name),
                        INDEX idx_proteins_created (created_at),
                        
                        FULLTEXT KEY ft_sequence (sequence),
                        FULLTEXT KEY ft_description (description)
                    ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci
                """))
                
                session.commit()
                self.stats['tables_created'] += 1
                self.logger.info("Created proteins table")
                
            except Exception as e:
                self.logger.error(f"Failed to create proteins table: {e}")
                self.stats['errors'] += 1
                raise
    
    def update_views(self) -> None:
        """Update database views to use the correct table structure."""
        self.logger.info("Updating database views")
        
        if self.dry_run:
            self.logger.info("DRY RUN: Would update database views")
            return
        
        with self.db_manager.get_session() as session:
            try:
                # Update protein_summary view
                session.execute(text("""
                    CREATE OR REPLACE VIEW protein_summary AS
                    SELECT 
                        p.isoform_id,
                        p.parent_protein_id,
                        p.name as protein_name,
                        p.organism,
                        p.sequence_length,
                        p.exon_count,
                        p.parent_tim_barrel_accession as tim_barrel_accession,
                        tbe.name as tim_barrel_name,
                        tbe.entry_type as tim_barrel_type,
                        CASE WHEN p.tim_barrel_location IS NOT NULL THEN 1 ELSE 0 END as has_tim_barrel,
                        p.created_at
                    FROM proteins p
                    JOIN interpro_proteins ip ON p.parent_protein_id = ip.uniprot_id 
                        AND p.parent_tim_barrel_accession = ip.tim_barrel_accession
                    JOIN tim_barrel_entries tbe ON ip.tim_barrel_accession = tbe.accession
                """))
                
                # Update collection_stats view
                session.execute(text("""
                    CREATE OR REPLACE VIEW collection_stats AS
                    SELECT 
                        (SELECT COUNT(*) FROM tim_barrel_entries) as tim_barrel_entries_count,
                        (SELECT COUNT(*) FROM tim_barrel_entries WHERE entry_type = 'pfam') as pfam_entries_count,
                        (SELECT COUNT(*) FROM tim_barrel_entries WHERE entry_type = 'interpro') as interpro_entries_count,
                        (SELECT COUNT(*) FROM interpro_proteins) as interpro_proteins_count,
                        (SELECT COUNT(DISTINCT uniprot_id) FROM interpro_proteins) as unique_proteins_count,
                        (SELECT COUNT(*) FROM proteins) as protein_isoforms_count,
                        (SELECT COUNT(*) FROM proteins WHERE tim_barrel_location IS NOT NULL) as tim_barrel_proteins_count,
                        (SELECT AVG(sequence_length) FROM proteins) as avg_sequence_length,
                        (SELECT AVG(exon_count) FROM proteins WHERE exon_count IS NOT NULL) as avg_exon_count
                """))
                
                session.commit()
                self.logger.info("Updated database views")
                
            except Exception as e:
                self.logger.error(f"Failed to update views: {e}")
                self.stats['errors'] += 1
    
    def run_cleanup(self) -> None:
        """Run the complete database cleanup."""
        self.logger.info("Starting database schema cleanup")
        
        try:
            # Check current tables
            existing_tables = self.check_current_tables()
            
            # Drop redundant tables
            self.drop_redundant_tables(existing_tables)
            
            # Create correct tables
            self.create_interpro_proteins_table()
            self.create_proteins_table()
            
            # Update views
            self.update_views()
            
            self.logger.info("Database cleanup completed successfully")
            
        except Exception as e:
            self.logger.error(f"Cleanup failed: {e}")
            raise
    
    def print_summary(self) -> None:
        """Print cleanup summary statistics."""
        print("\n" + "="*70)
        print("DATABASE CLEANUP SUMMARY")
        print("="*70)
        print(f"Redundant tables dropped:          {self.stats['tables_dropped']}")
        print(f"New tables created:                {self.stats['tables_created']}")
        print(f"Total errors:                      {self.stats['errors']}")
        print("="*70)
        
        if self.stats['errors'] > 0:
            print(f"⚠️  {self.stats['errors']} errors occurred during cleanup")
        else:
            print("✅ Cleanup completed successfully!")
        
        if self.dry_run:
            print("🔍 DRY RUN MODE - No data was actually modified")


def main():
    """Main script entry point."""
    parser = argparse.ArgumentParser(
        description="Clean up redundant database tables and ensure correct schema"
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
        
        # Create cleaner and run cleanup
        cleaner = DatabaseCleaner(dry_run=args.dry_run)
        cleaner.run_cleanup()
        
        # Print summary
        cleaner.print_summary()
        
    except KeyboardInterrupt:
        logger.info("Cleanup interrupted by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Cleanup failed: {str(e)}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()