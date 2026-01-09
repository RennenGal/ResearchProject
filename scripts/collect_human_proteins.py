#!/usr/bin/env python3
"""
Script to collect all Homo sapiens proteins for each TIM barrel entry found in step 1.

This script:
1. Reads all TIM barrel entries from the database (both PFAM families and InterPro entries)
2. For each PFAM family, queries InterPro API to get all human proteins
3. For each InterPro entry, queries InterPro API to get all human proteins
4. Stores the protein data in the database with proper relationships
5. Provides progress reporting and error handling

Usage:
    python scripts/collect_human_proteins.py [--config CONFIG_FILE] [--dry-run] [--verbose]
"""

import asyncio
import argparse
import logging
import sys
from pathlib import Path
from typing import List, Dict, Any, Set
from datetime import datetime

# Add the project root to Python path
sys.path.insert(0, str(Path(__file__).parent.parent))

from protein_data_collector.config import load_config_from_file, set_config, LoggingConfig
from protein_data_collector.api.interpro_client import InterProAPIClient
from protein_data_collector.models.entities import InterProProteinModel
from protein_data_collector.database.connection import get_database_manager
from protein_data_collector.database.schema import TIMBarrelEntry, InterProProtein
from protein_data_collector.logging_config import setup_logging


class HumanProteinCollector:
    """Collects human proteins for all TIM barrel entries from InterPro."""
    
    def __init__(self, dry_run: bool = False):
        """
        Initialize the collector.
        
        Args:
            dry_run: If True, don't actually write to database
        """
        self.dry_run = dry_run
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
        self.db_manager = get_database_manager()
        
        # Statistics
        self.stats = {
            'tim_barrel_entries_processed': 0,
            'pfam_entries_processed': 0,
            'interpro_entries_processed': 0,
            'proteins_found': 0,
            'proteins_stored': 0,
            'proteins_updated': 0,
            'proteins_skipped': 0,
            'errors': 0
        }
        
        # Track processed proteins to avoid duplicates
        self.processed_proteins: Set[str] = set()
    
    def get_tim_barrel_entries(self) -> List[TIMBarrelEntry]:
        """
        Get all TIM barrel entries from the database.
        
        Returns:
            List of TIMBarrelEntry instances from database
        """
        self.logger.info("Loading TIM barrel entries from database")
        
        with self.db_manager.get_session() as session:
            entries = session.query(TIMBarrelEntry).order_by(TIMBarrelEntry.accession).all()
            
            # Separate by type for logging
            pfam_entries = [e for e in entries if e.is_pfam]
            interpro_entries = [e for e in entries if e.is_interpro]
            
            self.logger.info(
                f"Loaded {len(entries)} TIM barrel entries: {len(pfam_entries)} PFAM + {len(interpro_entries)} InterPro",
                extra={
                    'total_entries': len(entries),
                    'pfam_count': len(pfam_entries),
                    'interpro_count': len(interpro_entries)
                }
            )
            
            return entries
    
    async def collect_proteins_for_pfam_entry(
        self, 
        client: InterProAPIClient, 
        entry: TIMBarrelEntry,
        page_size: int = 200
    ) -> List[InterProProteinModel]:
        """
        Collect human proteins for a PFAM family entry.
        
        Args:
            client: InterPro API client
            entry: PFAM family entry
            page_size: Number of results per API page
            
        Returns:
            List of validated InterProProteinModel instances
        """
        self.logger.info(f"Collecting proteins for PFAM family: {entry.accession} - {entry.name}")
        
        try:
            # Get proteins from InterPro API
            proteins_data = await client.get_proteins_in_pfam_family(
                entry.accession,
                organism="Homo sapiens",
                page_size=page_size
            )
            
            proteins = []
            for protein_data in proteins_data:
                try:
                    # Parse protein data
                    protein = client.parse_protein_data(protein_data, entry.accession)
                    
                    # Check for duplicates
                    if protein.uniprot_id not in self.processed_proteins:
                        proteins.append(protein)
                        self.processed_proteins.add(protein.uniprot_id)
                    else:
                        self.logger.debug(f"Skipping duplicate protein: {protein.uniprot_id}")
                        self.stats['proteins_skipped'] += 1
                        
                except Exception as e:
                    self.logger.warning(f"Failed to parse protein data for {entry.accession}: {str(e)}")
                    self.stats['errors'] += 1
            
            self.logger.info(
                f"Found {len(proteins)} unique human proteins for PFAM {entry.accession}",
                extra={
                    'pfam_accession': entry.accession,
                    'proteins_found': len(proteins),
                    'total_proteins_data': len(proteins_data)
                }
            )
            
            return proteins
            
        except Exception as e:
            self.logger.error(f"Failed to collect proteins for PFAM {entry.accession}: {str(e)}")
            self.stats['errors'] += 1
            return []
    
    async def collect_proteins_for_interpro_entry(
        self, 
        client: InterProAPIClient, 
        entry: TIMBarrelEntry,
        page_size: int = 200
    ) -> List[InterProProteinModel]:
        """
        Collect human proteins for an InterPro entry.
        
        Args:
            client: InterPro API client
            entry: InterPro entry
            page_size: Number of results per API page
            
        Returns:
            List of validated InterProProteinModel instances
        """
        self.logger.info(f"Collecting proteins for InterPro entry: {entry.accession} - {entry.name}")
        
        try:
            # Get proteins from InterPro API using the InterPro entry endpoint
            proteins_data = await client.get_proteins_in_interpro_entry(
                entry.accession,
                organism="Homo sapiens",
                page_size=page_size
            )
            
            proteins = []
            for protein_data in proteins_data:
                try:
                    # Parse protein data (adapt for InterPro entry)
                    protein = self._parse_interpro_protein_data(protein_data, entry.accession)
                    
                    # Check for duplicates
                    if protein.uniprot_id not in self.processed_proteins:
                        proteins.append(protein)
                        self.processed_proteins.add(protein.uniprot_id)
                    else:
                        self.logger.debug(f"Skipping duplicate protein: {protein.uniprot_id}")
                        self.stats['proteins_skipped'] += 1
                        
                except Exception as e:
                    self.logger.warning(f"Failed to parse protein data for {entry.accession}: {str(e)}")
                    self.stats['errors'] += 1
            
            self.logger.info(
                f"Found {len(proteins)} unique human proteins for InterPro {entry.accession}",
                extra={
                    'interpro_accession': entry.accession,
                    'proteins_found': len(proteins),
                    'total_proteins_data': len(proteins_data)
                }
            )
            
            return proteins
            
        except Exception as e:
            self.logger.error(f"Failed to collect proteins for InterPro {entry.accession}: {str(e)}")
            self.stats['errors'] += 1
            return []
    
    def _parse_interpro_protein_data(self, protein_data: Dict[str, Any], interpro_accession: str) -> InterProProteinModel:
        """
        Parse InterPro API response data into InterProProteinModel for InterPro entries.
        
        Args:
            protein_data: Raw protein data from InterPro API
            interpro_accession: Associated InterPro entry accession
            
        Returns:
            Validated InterProProteinModel instance
        """
        # Extract UniProt ID from metadata
        uniprot_id = protein_data.get('metadata', {}).get('accession')
        if not uniprot_id:
            raise ValueError("Missing UniProt ID in protein data")
        
        # Extract protein name and organism
        name = protein_data.get('metadata', {}).get('name', '')
        
        # Extract organism information
        organism = "Homo sapiens"  # Default, as we filter by this
        source_organism = protein_data.get('metadata', {}).get('source_organism', {})
        if isinstance(source_organism, dict):
            organism_name = source_organism.get('fullName', source_organism.get('scientificName', ''))
            if organism_name:
                organism = organism_name
        
        # Extract additional metadata
        basic_metadata = {
            'source_database': protein_data.get('metadata', {}).get('source_database', ''),
            'length': protein_data.get('metadata', {}).get('length'),
            'gene_name': protein_data.get('metadata', {}).get('gene', {}).get('name', ''),
            'protein_existence': protein_data.get('metadata', {}).get('protein_existence')
        }
        
        # Remove None values from metadata
        basic_metadata = {k: v for k, v in basic_metadata.items() if v is not None}
        
        return InterProProteinModel(
            uniprot_id=uniprot_id,
            tim_barrel_accession=interpro_accession,  # Use InterPro accession
            name=name,
            organism=organism,
            basic_metadata=basic_metadata
        )
    
    async def collect_all_human_proteins(self, page_size: int = 200) -> List[InterProProteinModel]:
        """
        Collect human proteins for all TIM barrel entries.
        
        Args:
            page_size: Number of results per API page
            
        Returns:
            List of all collected InterProProteinModel instances
        """
        self.logger.info("Starting comprehensive human protein collection for all TIM barrel entries")
        
        # Get all TIM barrel entries from database
        entries = self.get_tim_barrel_entries()
        
        if not entries:
            self.logger.warning("No TIM barrel entries found in database")
            return []
        
        all_proteins = []
        
        async with InterProAPIClient() as client:
            for entry in entries:
                try:
                    self.logger.info(f"Processing entry {entry.accession} ({entry.entry_type}): {entry.name}")
                    
                    proteins = []
                    if entry.is_pfam:
                        # Collect proteins for PFAM family
                        proteins = await self.collect_proteins_for_pfam_entry(client, entry, page_size)
                        self.stats['pfam_entries_processed'] += 1
                    elif entry.is_interpro:
                        # Collect proteins for InterPro entry
                        proteins = await self.collect_proteins_for_interpro_entry(client, entry, page_size)
                        self.stats['interpro_entries_processed'] += 1
                    
                    all_proteins.extend(proteins)
                    self.stats['proteins_found'] += len(proteins)
                    self.stats['tim_barrel_entries_processed'] += 1
                    
                    self.logger.info(
                        f"Completed processing {entry.accession}: found {len(proteins)} proteins",
                        extra={
                            'entry_accession': entry.accession,
                            'entry_type': entry.entry_type,
                            'proteins_found': len(proteins),
                            'total_proteins_so_far': len(all_proteins)
                        }
                    )
                    
                except Exception as e:
                    self.logger.error(f"Failed to process entry {entry.accession}: {str(e)}")
                    self.stats['errors'] += 1
        
        self.logger.info(
            f"Completed human protein collection: {len(all_proteins)} total proteins from {len(entries)} TIM barrel entries",
            extra={
                'total_proteins': len(all_proteins),
                'total_entries': len(entries),
                'pfam_entries': self.stats['pfam_entries_processed'],
                'interpro_entries': self.stats['interpro_entries_processed']
            }
        )
        
        return all_proteins
    
    def store_proteins(self, proteins: List[InterProProteinModel]) -> None:
        """
        Store human proteins in the database.
        
        Args:
            proteins: List of validated InterProProteinModel instances
        """
        if self.dry_run:
            self.logger.info(f"DRY RUN: Would store {len(proteins)} human proteins to database")
            self.stats['proteins_stored'] = len(proteins)
            return
        
        self.logger.info(f"Storing {len(proteins)} human proteins to database")
        
        try:
            with self.db_manager.get_transaction() as session:
                for protein in proteins:
                    try:
                        # Check if protein already exists
                        existing = session.query(InterProProtein).filter_by(
                            uniprot_id=protein.uniprot_id,
                            tim_barrel_accession=protein.tim_barrel_accession
                        ).first()
                        
                        if existing:
                            self.logger.debug(f"Updating existing protein: {protein.uniprot_id}")
                            # Update existing record
                            existing.name = protein.name
                            existing.organism = protein.organism
                            self.stats['proteins_updated'] += 1
                        else:
                            self.logger.debug(f"Creating new protein: {protein.uniprot_id}")
                            # Create new record
                            db_protein = InterProProtein(
                                uniprot_id=protein.uniprot_id,
                                tim_barrel_accession=protein.tim_barrel_accession,
                                name=protein.name,
                                organism=protein.organism,
                                created_at=datetime.now()
                            )
                            session.add(db_protein)
                            self.stats['proteins_stored'] += 1
                        
                    except Exception as e:
                        self.logger.error(
                            f"Failed to store protein {protein.uniprot_id}: {str(e)}",
                            extra={'protein': protein.model_dump()}
                        )
                        self.stats['errors'] += 1
                        # Continue with other proteins
                
                total_processed = self.stats['proteins_stored'] + self.stats['proteins_updated']
                self.logger.info(f"Successfully processed {total_processed} proteins ({self.stats['proteins_stored']} new, {self.stats['proteins_updated']} updated)")
                
        except Exception as e:
            self.logger.error(f"Database transaction failed: {str(e)}")
            raise
    
    def print_summary(self) -> None:
        """Print collection summary statistics."""
        print("\n" + "="*70)
        print("HUMAN PROTEIN COLLECTION SUMMARY")
        print("="*70)
        print(f"TIM barrel entries processed:      {self.stats['tim_barrel_entries_processed']}")
        print(f"  • PFAM entries processed:        {self.stats['pfam_entries_processed']}")
        print(f"  • InterPro entries processed:    {self.stats['interpro_entries_processed']}")
        print(f"Human proteins found:              {self.stats['proteins_found']}")
        print(f"Proteins skipped (duplicates):     {self.stats['proteins_skipped']}")
        print(f"New proteins stored:               {self.stats['proteins_stored']}")
        print(f"Existing proteins updated:         {self.stats['proteins_updated']}")
        print(f"Total errors:                      {self.stats['errors']}")
        print("="*70)
        
        if self.stats['errors'] > 0:
            print(f"⚠️  {self.stats['errors']} errors occurred during collection")
        else:
            print("✅ Collection completed successfully!")
        
        if self.dry_run:
            print("🔍 DRY RUN MODE - No data was actually stored")


async def main():
    """Main script entry point."""
    parser = argparse.ArgumentParser(
        description="Collect all Homo sapiens proteins for TIM barrel entries from InterPro"
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
        help='Run without storing data to database'
    )
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose logging'
    )
    parser.add_argument(
        '--page-size',
        type=int,
        default=200,
        help='Number of results per API page (default: 200)'
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
        
        # Create collector and run collection
        collector = HumanProteinCollector(dry_run=args.dry_run)
        
        # Collect all human proteins
        proteins = await collector.collect_all_human_proteins(page_size=args.page_size)
        
        if not proteins:
            logger.warning("No human proteins found")
            return
        
        # Store proteins to database
        collector.store_proteins(proteins)
        
        # Print summary
        collector.print_summary()
        
        # Show some example proteins
        if args.verbose and proteins:
            print(f"\nExample proteins collected ({len(proteins)} total):")
            for i, protein in enumerate(proteins[:5], 1):
                print(f"\n{i}. {protein.uniprot_id} - {protein.name}")
                print(f"   TIM barrel entry: {protein.tim_barrel_accession}")
                print(f"   Organism: {protein.organism}")
                if protein.basic_metadata:
                    print(f"   Metadata: {protein.basic_metadata}")
        
    except KeyboardInterrupt:
        logger.info("Collection interrupted by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Collection failed: {str(e)}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    asyncio.run(main())