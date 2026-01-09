#!/usr/bin/env python3
"""
Summary of TIM barrel PFAM family collection status.
"""

import sys
from pathlib import Path

# Add the project root to Python path
sys.path.insert(0, str(Path(__file__).parent.parent))

from protein_data_collector.config import load_config_from_file, set_config
from protein_data_collector.database.connection import get_database_manager
from protein_data_collector.database.schema import TIMBarrelEntry

def show_tim_barrel_summary():
    """Show current status of TIM barrel entry collection from unified table."""
    
    # Load configuration
    config_path = Path("config.test.json")
    config = load_config_from_file(str(config_path))
    set_config(config)
    
    # Get database manager
    db_manager = get_database_manager()
    
    print("🧬 TIM BARREL COLLECTION STATUS")
    print("=" * 70)
    
    with db_manager.get_session() as session:
        # Get all entries from unified table
        all_entries = session.query(TIMBarrelEntry).order_by(TIMBarrelEntry.accession).all()
        
        # Separate by type
        pfam_entries = [entry for entry in all_entries if entry.is_pfam]
        interpro_entries = [entry for entry in all_entries if entry.is_interpro]
        
        print(f"📊 CURRENT STATUS:")
        print(f"   • PFAM families collected: {len(pfam_entries)}")
        print(f"   • InterPro entries collected: {len(interpro_entries)}")
        print(f"   • Total TIM barrel entries: {len(all_entries)}")
        print(f"   • Using unified table storage")
        print(f"   • Comprehensive hybrid search strategy")
        print(f"   • No artificial limits - collecting all TIM barrel entries found")
        print()
        
        if pfam_entries:
            print(f"✅ COLLECTED PFAM FAMILIES ({len(pfam_entries)}):")
            for i, entry in enumerate(pfam_entries, 1):
                print(f"   {i:2d}. {entry.accession} - {entry.name}")
                if entry.interpro_id:
                    print(f"       InterPro: {entry.interpro_id}")
            print()
        
        if interpro_entries:
            print(f"✅ COLLECTED INTERPRO ENTRIES ({len(interpro_entries)}):")
            for i, entry in enumerate(interpro_entries, 1):
                print(f"   {i:2d}. {entry.accession} - {entry.name}")
                if entry.interpro_type:
                    print(f"       Type: {entry.interpro_type}")
            print()
        
        print("🔍 COMPREHENSIVE HYBRID SEARCH STRATEGY:")
        print("   • PHASE 1: Direct PFAM family search with multiple terms")
        print("   • PHASE 2: InterPro entry (IPR) search for structural classifications")
        print("   • Automatic deduplication of results")
        print("   • Unified storage in single table")
        print("   • No artificial limits on entry count")
        print()
        
        print("💡 SEARCH COVERAGE:")
        print("   • Text-based searches for TIM barrel variants")
        print("   • Structural classification searches (aldolase-type, etc.)")
        print("   • InterPro domain and family classifications")
        print("   • Known TIM barrel enzyme families")
        print()
        
        if len(all_entries) >= 49:
            print("🎉 COMPREHENSIVE COLLECTION COMPLETE!")
            print(f"   We have collected {len(all_entries)} TIM barrel entries")
            print("   Our hybrid search strategy successfully found both PFAM families and InterPro entries")
        elif len(all_entries) >= 22:
            print("🎉 EXCEEDED INITIAL TARGET!")
            print(f"   We have {len(all_entries)} entries, exceeding the initial 22 family target")
            print("   Our hybrid search strategy is finding comprehensive TIM barrel classifications")
        else:
            print("📈 COLLECTION IN PROGRESS")
            print("   Run the unified collection script to gather all available TIM barrel entries")
            print("   python scripts/collect_tim_barrel_entries.py")

if __name__ == "__main__":
    show_tim_barrel_summary()