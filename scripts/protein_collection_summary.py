#!/usr/bin/env python3
"""
Summary of human protein collection status for TIM barrel entries.
"""

import sys
from pathlib import Path

# Add the project root to Python path
sys.path.insert(0, str(Path(__file__).parent.parent))

from protein_data_collector.config import load_config_from_file, set_config
from protein_data_collector.database.connection import get_database_manager
from protein_data_collector.database.schema import TIMBarrelEntry, InterProProtein

def show_protein_collection_summary():
    """Show current status of human protein collection for TIM barrel entries."""
    
    # Load configuration
    config_path = Path("config.test.json")
    config = load_config_from_file(str(config_path))
    set_config(config)
    
    # Get database manager
    db_manager = get_database_manager()
    
    print("🧬 HUMAN PROTEIN COLLECTION STATUS")
    print("=" * 70)
    
    with db_manager.get_session() as session:
        # Get all TIM barrel entries
        all_entries = session.query(TIMBarrelEntry).order_by(TIMBarrelEntry.accession).all()
        
        # Get all collected proteins
        all_proteins = session.query(InterProProtein).all()
        
        # Separate entries by type
        pfam_entries = [entry for entry in all_entries if entry.is_pfam]
        interpro_entries = [entry for entry in all_entries if entry.is_interpro]
        
        # Group proteins by TIM barrel entry
        proteins_by_entry = {}
        for protein in all_proteins:
            if protein.tim_barrel_accession not in proteins_by_entry:
                proteins_by_entry[protein.tim_barrel_accession] = []
            proteins_by_entry[protein.tim_barrel_accession].append(protein)
        
        print(f"📊 OVERALL STATUS:")
        print(f"   • TIM barrel entries: {len(all_entries)} (18 PFAM + 31 InterPro)")
        print(f"   • Entries with proteins collected: {len(proteins_by_entry)}")
        print(f"   • Total human proteins collected: {len(all_proteins)}")
        print(f"   • Unique proteins: {len(set(p.uniprot_id for p in all_proteins))}")
        print()
        
        # Show collection progress
        entries_with_proteins = len(proteins_by_entry)
        collection_progress = (entries_with_proteins / len(all_entries)) * 100 if all_entries else 0
        
        print(f"📈 COLLECTION PROGRESS:")
        print(f"   • Progress: {collection_progress:.1f}% ({entries_with_proteins}/{len(all_entries)} entries)")
        
        if entries_with_proteins == 0:
            print("   • Status: Not started - run collect_human_proteins.py")
        elif entries_with_proteins < len(all_entries):
            print("   • Status: In progress - some entries still need processing")
        else:
            print("   • Status: Complete - all entries processed")
        print()
        
        # Show detailed breakdown
        if proteins_by_entry:
            print(f"🔍 DETAILED BREAKDOWN:")
            
            # PFAM entries with proteins
            pfam_with_proteins = [e for e in pfam_entries if e.accession in proteins_by_entry]
            if pfam_with_proteins:
                print(f"\n✅ PFAM FAMILIES WITH PROTEINS ({len(pfam_with_proteins)}):")
                for entry in pfam_with_proteins:
                    protein_count = len(proteins_by_entry[entry.accession])
                    print(f"   {entry.accession} - {entry.name}: {protein_count} proteins")
            
            # InterPro entries with proteins
            interpro_with_proteins = [e for e in interpro_entries if e.accession in proteins_by_entry]
            if interpro_with_proteins:
                print(f"\n✅ INTERPRO ENTRIES WITH PROTEINS ({len(interpro_with_proteins)}):")
                for entry in interpro_with_proteins:
                    protein_count = len(proteins_by_entry[entry.accession])
                    print(f"   {entry.accession} - {entry.name}: {protein_count} proteins")
            
            # Show entries without proteins
            entries_without_proteins = [e for e in all_entries if e.accession not in proteins_by_entry]
            if entries_without_proteins:
                print(f"\n⏳ ENTRIES PENDING PROTEIN COLLECTION ({len(entries_without_proteins)}):")
                for entry in entries_without_proteins[:10]:  # Show first 10
                    print(f"   {entry.accession} ({entry.entry_type}) - {entry.name}")
                if len(entries_without_proteins) > 10:
                    print(f"   ... and {len(entries_without_proteins) - 10} more")
        
        # Show top protein-rich entries
        if proteins_by_entry:
            print(f"\n🏆 TOP ENTRIES BY PROTEIN COUNT:")
            sorted_entries = sorted(proteins_by_entry.items(), key=lambda x: len(x[1]), reverse=True)
            for i, (accession, proteins) in enumerate(sorted_entries[:5], 1):
                entry = session.query(TIMBarrelEntry).filter_by(accession=accession).first()
                entry_name = entry.name if entry else "Unknown"
                print(f"   {i}. {accession} - {entry_name}: {len(proteins)} proteins")
        
        # Show some example proteins
        if all_proteins:
            print(f"\n🧪 EXAMPLE COLLECTED PROTEINS:")
            for i, protein in enumerate(all_proteins[:5], 1):
                print(f"   {i}. {protein.uniprot_id} - {protein.name or 'No name'}")
                print(f"      TIM barrel entry: {protein.tim_barrel_accession}")
                print(f"      Organism: {protein.organism}")
        
        print()
        
        # Show next steps
        if entries_with_proteins == 0:
            print("🚀 NEXT STEPS:")
            print("   Run: python scripts/collect_human_proteins.py")
            print("   This will collect all human proteins for each TIM barrel entry")
        elif entries_with_proteins < len(all_entries):
            print("🚀 NEXT STEPS:")
            print("   Continue running: python scripts/collect_human_proteins.py")
            print("   Or check logs if collection is in progress")
        else:
            print("🎉 COLLECTION COMPLETE!")
            print("   All TIM barrel entries have been processed")
            print("   Ready for next analysis steps")

if __name__ == "__main__":
    show_protein_collection_summary()