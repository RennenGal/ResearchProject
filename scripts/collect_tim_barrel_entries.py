#!/usr/bin/env python3
"""
Unified script to collect both PFAM families and InterPro entries with TIM barrel annotations
from InterPro API and upload them to the database.

This script:
1. Queries InterPro API for both PFAM families and InterPro entries with TIM barrel annotations
2. Parses and validates the response data
3. Stores all entries in a single unified table
4. Provides progress reporting and error handling

Usage:
    python scripts/collect_tim_barrel_entries.py [--config CONFIG_FILE] [--dry-run] [--verbose]
"""

import asyncio
import argparse
import logging
import sys
from pathlib import Path
from typing import List, Dict, Any
from datetime import datetime

# Add the project root to Python path
sys.path.insert(0, str(Path(__file__).parent.parent))

from protein_data_collector.config import load_config_from_file, set_config, LoggingConfig
from protein_data_collector.api.interpro_client import InterProAPIClient
from protein_data_collector.models.entities import TIMBarrelEntryModel
from protein_data_collector.database.connection import get_database_manager
from protein_data_collector.database.schema import TIMBarrelEntry
from protein_data_collector.logging_config import setup_logging


class UnifiedTIMBarrelCollector:
    """Collects both PFAM families and InterPro entries with TIM barrel annotations from InterPro."""
    
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
            'pfam_entries_found': 0,
            'interpro_entries_found': 0,
            'total_entries_found': 0,
            'entries_stored': 0,
            'entries_updated': 0,
            'errors': 0
        }
    
    async def collect_tim_barrel_entries(self, page_size: int = 200) -> List[TIMBarrelEntryModel]:
        """
        Collect both PFAM families and InterPro entries with TIM barrel annotations from InterPro.
        
        This method uses a comprehensive approach:
        1. Direct PFAM family search with multiple search terms
        2. InterPro entry (IPR) search to find structural classifications
        3. Deduplicate results and store in unified format
        
        Args:
            page_size: Number of results per API page
            
        Returns:
            List of TIMBarrelEntryModel instances
        """
        self.logger.info("Starting comprehensive TIM barrel collection from InterPro API")
        
        all_entries = {}  # Use dict to avoid duplicates, keyed by accession
        
        async with InterProAPIClient() as client:
            # PHASE 1: Direct PFAM family search
            await self._collect_pfam_families(client, all_entries, page_size)
            
            # PHASE 2: InterPro entry search
            await self._collect_interpro_entries(client, all_entries, page_size)
        
        # Convert to list and update stats
        entries_list = list(all_entries.values())
        
        # Count by type
        pfam_count = sum(1 for entry in entries_list if entry.is_pfam)
        interpro_count = sum(1 for entry in entries_list if entry.is_interpro)
        
        self.stats['pfam_entries_found'] = pfam_count
        self.stats['interpro_entries_found'] = interpro_count
        self.stats['total_entries_found'] = len(entries_list)
        
        self.logger.info(f"Completed comprehensive collection: {pfam_count} PFAM families + {interpro_count} InterPro entries = {len(entries_list)} total TIM barrel entries")
        return entries_list
    
    async def _collect_pfam_families(self, client: InterProAPIClient, all_entries: dict, page_size: int):
        """Collect PFAM families using direct search terms."""
        self.logger.info("PHASE 1: Direct PFAM family search")
        
        # Use multiple search terms to catch all TIM barrel variants
        search_terms = [
            "TIM barrel",
            "TIM-barrel", 
            "tim-barrel",
            "TIM",  # Catches TIM-barrel, tim-barrel, etc.
            "barrel",  # Catches other barrel variants
            "triosephosphate",  # Catches the classic TIM barrel protein
            "aldolase",  # Aldolase-type TIM barrels
            "isomerase",  # Many TIM barrel enzymes are isomerases
            "glycosyl hydrolase",  # Several TIM barrel glycosyl hydrolases
            "enolase",  # Enolase has TIM barrel domain
            "synthase",  # Some synthases have TIM barrel domains
        ]
        
        for search_term in search_terms:
            try:
                self.logger.info(f"Searching PFAM families with term: '{search_term}'")
                
                # Get raw family data from InterPro
                families_data = await client.get_pfam_families_with_tim_barrel_search(search_term, page_size)
                
                self.logger.info(f"Found {len(families_data)} PFAM families for search term '{search_term}'")
                
                # Filter and collect TIM barrel families
                for family_data in families_data:
                    try:
                        metadata = family_data.get('metadata', {})
                        accession = metadata.get('accession', '')
                        name = metadata.get('name', '').lower()
                        
                        # Check if it's actually TIM barrel related
                        is_tim_barrel = self._is_tim_barrel_family(name, accession)
                        
                        if is_tim_barrel and accession:
                            # Parse the family data into a validated model
                            entry = client.parse_pfam_family_data(family_data)
                            all_entries[accession] = entry
                            
                            self.logger.debug(
                                f"Found TIM barrel PFAM family: {entry.accession} - {entry.name}",
                                extra={
                                    'accession': entry.accession,
                                    'name': entry.name,
                                    'search_term': search_term,
                                    'source': 'direct_pfam_search'
                                }
                            )
                    
                    except Exception as e:
                        self.logger.warning(f"Failed to process PFAM family data: {str(e)}")
                        self.stats['errors'] += 1
            
            except Exception as e:
                self.logger.error(f"Failed to search PFAM families with term '{search_term}': {str(e)}")
                self.stats['errors'] += 1
        
        pfam_count = sum(1 for entry in all_entries.values() if entry.is_pfam)
        self.logger.info(f"PHASE 1 completed: {pfam_count} unique PFAM families from direct search")
    
    async def _collect_interpro_entries(self, client: InterProAPIClient, all_entries: dict, page_size: int):
        """Collect InterPro entries (IPR records) with TIM barrel annotations."""
        self.logger.info("PHASE 2: InterPro entry collection")
        
        # Search terms for InterPro entries
        interpro_search_terms = [
            "TIM barrel",
            "TIM-barrel",
            "aldolase TIM barrel",
            "triosephosphate isomerase",
            "barrel fold",
            "eight-stranded beta/alpha barrel",
            "8-stranded beta/alpha barrel"
        ]
        
        # Search for InterPro entries
        for search_term in interpro_search_terms:
            try:
                self.logger.info(f"Searching InterPro entries with term: '{search_term}'")
                
                entries_data = await client.get_interpro_entries_with_tim_barrel_search(search_term, page_size)
                
                self.logger.info(f"Found {len(entries_data)} InterPro entries for search term '{search_term}'")
                
                # Filter and collect TIM barrel InterPro entries
                for entry_data in entries_data:
                    try:
                        metadata = entry_data.get('metadata', {})
                        accession = metadata.get('accession', '')
                        name = metadata.get('name', '').lower()
                        
                        # Check if it's TIM barrel related
                        is_tim_barrel_entry = (
                            'tim' in name or 'barrel' in name or 'aldolase' in name or
                            'triosephosphate' in name or 'isomerase' in name or
                            accession.startswith('IPR013785') or  # Aldolase-type TIM barrel
                            'fold' in name or 'glycosyl' in name or 'hydrolase' in name
                        )
                        
                        if is_tim_barrel_entry and accession:
                            # Parse the entry data into a validated model
                            entry = client.parse_interpro_entry_data(entry_data)
                            all_entries[accession] = entry
                            
                            self.logger.debug(
                                f"Found TIM barrel InterPro entry: {entry.accession} - {entry.name}",
                                extra={
                                    'accession': entry.accession,
                                    'name': entry.name,
                                    'search_term': search_term,
                                    'source': 'interpro_entry_search'
                                }
                            )
                    
                    except Exception as e:
                        self.logger.warning(f"Failed to process InterPro entry data: {str(e)}")
                        self.stats['errors'] += 1
            
            except Exception as e:
                self.logger.error(f"Failed to search InterPro entries with term '{search_term}': {str(e)}")
                self.stats['errors'] += 1
        
        interpro_count = sum(1 for entry in all_entries.values() if entry.is_interpro)
        self.logger.info(f"PHASE 2 completed: {interpro_count} unique InterPro entries collected")
    
    def _is_tim_barrel_family(self, name: str, accession: str) -> bool:
        """Check if a family is TIM barrel related based on name and accession."""
        return (
            # Direct TIM barrel mentions
            ('tim' in name and 'barrel' in name) or
            ('tim-barrel' in name) or
            ('tim_barrel' in name) or
            ('tim barrel' in name) or
            
            # Classic triosephosphate isomerase (the original TIM barrel)
            ('triosephosphate' in name and 'isomerase' in name) or
            (accession == 'PF00121') or  # Explicit inclusion of classic TIM barrel
            
            # Barrel with TIM-related terms
            ('barrel' in name and any(word in name for word in ['tim', 'triosephosphate', 'isomerase'])) or
            
            # TIM with structural terms (but exclude non-barrel TIM proteins)
            ('tim' in name and 'isomerase' in name) or
            ('tim' in name and 'fold' in name and 'barrel' in name) or
            
            # Aldolase-type TIM barrels (common structural family)
            ('aldolase' in name and any(word in name for word in ['tim', 'barrel', 'fold'])) or
            
            # Specific enzyme families known to have TIM barrel structure
            ('enolase' in name and 'tim' in name) or
            ('xylose isomerase' in name) or
            ('malate synthase' in name and 'tim' in name) or
            ('glycosyl hydrolase' in name and 'tim' in name) or
            ('lactonase' in name and 'tim' in name) or
            ('phosphorylase' in name and 'tim' in name) or
            ('aminomutase' in name and 'tim' in name) or
            ('lyase' in name and 'tim' in name) or
            ('acetylglucosaminidase' in name and 'tim' in name) or
            ('prmt5' in name and 'tim' in name) or
            
            # Additional patterns that might be TIM barrels
            ('gcpe' in name and 'tim' in name) or
            ('gta' in name and 'tim' in name) or
            ('mtc6' in name and 'tim' in name) or
            ('endos' in name and 'tim' in name)
        )
    
    def store_entries(self, entries: List[TIMBarrelEntryModel]) -> None:
        """
        Store TIM barrel entries in the unified database table.
        
        Args:
            entries: List of validated TIMBarrelEntryModel instances
        """
        if self.dry_run:
            self.logger.info(f"DRY RUN: Would store {len(entries)} TIM barrel entries to database")
            self.stats['entries_stored'] = len(entries)
            return
        
        self.logger.info(f"Storing {len(entries)} TIM barrel entries to database")
        
        try:
            with self.db_manager.get_transaction() as session:
                for entry in entries:
                    try:
                        # Check if entry already exists
                        existing = session.query(TIMBarrelEntry).filter_by(
                            accession=entry.accession
                        ).first()
                        
                        if existing:
                            self.logger.debug(f"Updating existing entry: {entry.accession}")
                            # Update existing record
                            existing.entry_type = entry.entry_type
                            existing.name = entry.name
                            existing.description = entry.description
                            existing.interpro_type = entry.interpro_type
                            existing.tim_barrel_annotation = entry.tim_barrel_annotation
                            existing.member_databases = entry.member_databases
                            existing.interpro_id = entry.interpro_id
                            self.stats['entries_updated'] += 1
                        else:
                            self.logger.debug(f"Creating new entry: {entry.accession}")
                            # Create new record
                            db_entry = TIMBarrelEntry(
                                accession=entry.accession,
                                entry_type=entry.entry_type,
                                name=entry.name,
                                description=entry.description,
                                interpro_type=entry.interpro_type,
                                tim_barrel_annotation=entry.tim_barrel_annotation,
                                member_databases=entry.member_databases,
                                interpro_id=entry.interpro_id,
                                created_at=datetime.now()
                            )
                            session.add(db_entry)
                            self.stats['entries_stored'] += 1
                        
                    except Exception as e:
                        self.logger.error(
                            f"Failed to store entry {entry.accession}: {str(e)}",
                            extra={'entry': entry.model_dump()}
                        )
                        self.stats['errors'] += 1
                        # Continue with other entries
                
                total_processed = self.stats['entries_stored'] + self.stats['entries_updated']
                self.logger.info(f"Successfully processed {total_processed} entries ({self.stats['entries_stored']} new, {self.stats['entries_updated']} updated)")
                
        except Exception as e:
            self.logger.error(f"Database transaction failed: {str(e)}")
            raise
    
    def print_summary(self) -> None:
        """Print collection summary statistics."""
        print("\n" + "="*60)
        print("UNIFIED TIM BARREL COLLECTION SUMMARY")
        print("="*60)
        print(f"PFAM families found:           {self.stats['pfam_entries_found']}")
        print(f"InterPro entries found:        {self.stats['interpro_entries_found']}")
        print(f"Total entries found:           {self.stats['total_entries_found']}")
        print(f"New entries stored:            {self.stats['entries_stored']}")
        print(f"Existing entries updated:      {self.stats['entries_updated']}")
        print(f"Total errors:                  {self.stats['errors']}")
        print("="*60)
        
        if self.stats['errors'] > 0:
            print(f"⚠️  {self.stats['errors']} errors occurred during collection")
        else:
            print("✅ Collection completed successfully!")
        
        if self.dry_run:
            print("🔍 DRY RUN MODE - No data was actually stored")


async def main():
    """Main script entry point."""
    parser = argparse.ArgumentParser(
        description="Collect both PFAM families and InterPro entries with TIM barrel annotations from InterPro"
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
        collector = UnifiedTIMBarrelCollector(dry_run=args.dry_run)
        
        # Collect all TIM barrel entries
        entries = await collector.collect_tim_barrel_entries(page_size=args.page_size)
        
        if not entries:
            logger.warning("No TIM barrel entries found")
            return
        
        # Store entries to database
        collector.store_entries(entries)
        
        # Print summary
        collector.print_summary()
        
        # Show some example entries
        if args.verbose:
            pfam_entries = [e for e in entries if e.is_pfam]
            interpro_entries = [e for e in entries if e.is_interpro]
            
            if pfam_entries:
                print(f"\nExample PFAM families collected ({len(pfam_entries)} total):")
                for i, entry in enumerate(pfam_entries[:3], 1):
                    print(f"\n{i}. {entry.accession} - {entry.name}")
                    print(f"   Description: {entry.description[:100]}..." if entry.description else "   No description")
                    print(f"   TIM annotation: {entry.tim_barrel_annotation[:100]}...")
            
            if interpro_entries:
                print(f"\nExample InterPro entries collected ({len(interpro_entries)} total):")
                for i, entry in enumerate(interpro_entries[:3], 1):
                    print(f"\n{i}. {entry.accession} - {entry.name}")
                    print(f"   Type: {entry.interpro_type}")
                    print(f"   Description: {entry.description[:100]}..." if entry.description else "   No description")
                    print(f"   TIM annotation: {entry.tim_barrel_annotation[:100]}...")
        
    except KeyboardInterrupt:
        logger.info("Collection interrupted by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Collection failed: {str(e)}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    asyncio.run(main())