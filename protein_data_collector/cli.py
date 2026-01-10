"""
Command-line interface for the Protein Data Collector system.

This module provides CLI commands for data collection, querying, and export.
"""

import click
import sys
import requests
from pathlib import Path
from .config import SystemConfig, load_config_from_file
from .logging_config import setup_logging


@click.group()
@click.option('--config', '-c', type=click.Path(exists=True), 
              help='Configuration file path')
@click.option('--verbose', '-v', is_flag=True, 
              help='Enable verbose logging')
@click.pass_context
def cli(ctx, config, verbose):
    """Protein Data Collector - Bioinformatics data collection tool."""
    # Ensure context object exists
    ctx.ensure_object(dict)
    
    # Load configuration
    if config:
        system_config = load_config_from_file(config)
    else:
        system_config = SystemConfig.from_env()
    
    # Adjust logging level if verbose
    if verbose:
        system_config.logging.level = "DEBUG"
    
    # Setup logging
    setup_logging(system_config.logging)
    
    # Store config in context
    ctx.obj['config'] = system_config


@cli.command()
@click.option('--tim-barrel-query', default='TIM barrel', 
              help='Query string for TIM barrel annotations')
@click.option('--resume', is_flag=True, 
              help='Resume interrupted collection')
@click.option('--progress-file', type=click.Path(), 
              help='Progress file for resume functionality')
@click.option('--page-size', default=200, type=int,
              help='Number of results per page for API queries')
@click.option('--batch-size', type=int,
              help='Batch size for processing operations')
@click.option('--no-store', is_flag=True,
              help='Skip storing data in database (collect only)')
@click.pass_context
def collect(ctx, tim_barrel_query, resume, progress_file, page_size, batch_size, no_store):
    """Run the full data collection process."""
    import asyncio
    from .collector.data_collector import DataCollector, run_complete_collection, resume_collection
    
    config = ctx.obj['config']
    
    click.echo("=== Protein Data Collection ===")
    click.echo(f"TIM barrel query: {tim_barrel_query}")
    click.echo(f"Resume mode: {resume}")
    click.echo(f"Page size: {page_size}")
    click.echo(f"Batch size: {batch_size or config.collection.batch_size}")
    click.echo(f"Store data: {not no_store}")
    
    if progress_file:
        click.echo(f"Progress file: {progress_file}")
    
    try:
        if resume and progress_file:
            click.echo("\nResuming collection from progress file...")
            report = asyncio.run(resume_collection(progress_file))
        else:
            click.echo("\nStarting fresh collection...")
            report = asyncio.run(run_complete_collection(
                progress_file=progress_file,
                page_size=page_size,
                batch_size=batch_size,
                store_data=not no_store
            ))
        
        # Display results
        click.echo("\n=== Collection Results ===")
        click.echo(f"Phase completed: {report.progress.phase}")
        click.echo(f"Duration: {report.progress.duration_seconds:.1f} seconds")
        click.echo(f"PFAM families: {len(report.pfam_families)}")
        click.echo(f"InterPro proteins: {len(report.interpro_proteins)}")
        click.echo(f"UniProt isoforms: {len(report.uniprot_isoforms)}")
        click.echo(f"Success rate: {report.success_rate:.1f}%")
        
        if report.validation_errors or report.api_errors or report.storage_errors:
            click.echo("\n=== Errors ===")
            if report.validation_errors:
                click.echo(f"Validation errors: {len(report.validation_errors)}")
            if report.api_errors:
                click.echo(f"API errors: {len(report.api_errors)}")
            if report.storage_errors:
                click.echo(f"Storage errors: {len(report.storage_errors)}")
        
        if not no_store:
            click.echo(f"\nEntities stored: {report.progress.total_entities_stored}")
        
        click.echo("\nCollection completed successfully!")
        
    except KeyboardInterrupt:
        click.echo("\nCollection cancelled by user.")
        sys.exit(1)
    except Exception as e:
        click.echo(f"\nCollection failed: {e}", err=True)
        sys.exit(1)


@cli.command()
@click.option('--pfam-family', help='Search by PFAM family identifier')
@click.option('--protein-id', help='Search by protein identifier')
@click.option('--tim-barrel', is_flag=True, help='Search for proteins with TIM barrel annotations')
@click.option('--min-length', type=int, help='Minimum sequence length')
@click.option('--max-length', type=int, help='Maximum sequence length')
@click.option('--min-exons', type=int, help='Minimum exon count')
@click.option('--max-exons', type=int, help='Maximum exon count')
@click.option('--organism', default='Homo sapiens', help='Organism filter')
@click.option('--limit', default=100, type=int, help='Maximum number of results')
@click.option('--include-isoforms', is_flag=True, default=True, help='Include isoform details')
@click.option('--format', 'output_format', type=click.Choice(['table', 'json', 'summary']), 
              default='table', help='Output format')
@click.pass_context
def query(ctx, pfam_family, protein_id, tim_barrel, min_length, max_length, 
          min_exons, max_exons, organism, limit, include_isoforms, output_format):
    """Query collected protein data."""
    from .query.engine import QueryEngine, QueryFilters
    import json
    
    config = ctx.obj['config']
    engine = QueryEngine()
    
    click.echo("=== Protein Data Query ===")
    
    try:
        result = None
        
        if pfam_family:
            click.echo(f"Searching for PFAM family: {pfam_family}")
            result = engine.search_by_pfam_family(pfam_family, include_isoforms)
            
        elif protein_id:
            click.echo(f"Searching for protein: {protein_id}")
            result = engine.get_protein_isoforms(protein_id)
            
        elif tim_barrel:
            click.echo("Searching for proteins with TIM barrel annotations")
            criteria = {'has_location': True}
            result = engine.search_by_tim_barrel_features(criteria)
            
        else:
            # General filtering
            click.echo("Filtering proteins with specified criteria")
            filters = QueryFilters(
                organism=organism,
                min_sequence_length=min_length,
                max_sequence_length=max_length,
                min_exon_count=min_exons,
                max_exon_count=max_exons,
                has_tim_barrel=tim_barrel if tim_barrel else None
            )
            result = engine.filter_proteins(filters, limit)
        
        if not result or result.total_count == 0:
            click.echo("No results found.")
            return
        
        # Display results based on format
        if output_format == 'json':
            output_data = {
                'pfam_families': result.pfam_families,
                'proteins': result.proteins,
                'isoforms': result.isoforms if include_isoforms else [],
                'total_count': result.total_count,
                'metadata': result.query_metadata
            }
            click.echo(json.dumps(output_data, indent=2))
            
        elif output_format == 'summary':
            click.echo(f"\n=== Query Results Summary ===")
            click.echo(f"Total results: {result.total_count}")
            click.echo(f"PFAM families: {len(result.pfam_families)}")
            click.echo(f"Proteins: {len(result.proteins)}")
            if include_isoforms:
                click.echo(f"Isoforms: {len(result.isoforms)}")
            
        else:  # table format
            click.echo(f"\n=== Query Results ({result.total_count} total) ===")
            
            if result.pfam_families:
                click.echo("\nPFAM Families:")
                for family in result.pfam_families:
                    click.echo(f"  {family['accession']}: {family['name']}")
            
            if result.proteins:
                click.echo("\nProteins:")
                for protein in result.proteins[:limit]:
                    click.echo(f"  {protein['uniprot_id']}: {protein.get('name', 'N/A')}")
            
            if include_isoforms and result.isoforms:
                click.echo(f"\nIsoforms (showing first {min(10, len(result.isoforms))}):")
                for isoform in result.isoforms[:10]:
                    tim_info = ""
                    if isoform.get('tim_barrel_location'):
                        tim_loc = isoform['tim_barrel_location']
                        if isinstance(tim_loc, dict) and 'start' in tim_loc:
                            tim_info = f" [TIM: {tim_loc['start']}-{tim_loc['end']}]"
                    click.echo(f"  {isoform['isoform_id']}: {isoform['sequence_length']} aa{tim_info}")
        
    except Exception as e:
        click.echo(f"Query failed: {e}", err=True)
        sys.exit(1)


@cli.command()
@click.option('--format', type=click.Choice(['fasta', 'json', 'csv']), 
              default='json', help='Export format')
@click.option('--output', '-o', type=click.Path(), 
              help='Output file path (default: stdout)')
@click.option('--pfam-family', help='Filter by PFAM family')
@click.option('--protein-id', help='Filter by specific protein ID')
@click.option('--tim-barrel', is_flag=True, help='Filter proteins with TIM barrel annotations')
@click.option('--min-length', type=int, help='Minimum sequence length')
@click.option('--max-length', type=int, help='Maximum sequence length')
@click.option('--limit', type=int, help='Maximum number of results to export')
@click.option('--include-sequences/--no-sequences', default=True, 
              help='Include protein sequences in export')
@click.option('--include-metadata/--no-metadata', default=True,
              help='Include metadata (names, descriptions, etc.)')
@click.option('--include-annotations/--no-annotations', default=True,
              help='Include annotations (TIM barrel, exons, etc.)')
@click.option('--max-seq-length', type=int, help='Truncate sequences to this length')
@click.pass_context
def export(ctx, format, output, pfam_family, protein_id, tim_barrel, min_length, max_length, 
           limit, include_sequences, include_metadata, include_annotations, max_seq_length):
    """Export collected data in various formats."""
    from .query.engine import QueryEngine, QueryFilters
    from .query.export import DataExporter, ExportFormat, ExportOptions
    
    config = ctx.obj['config']
    engine = QueryEngine()
    exporter = DataExporter()
    
    click.echo(f"=== Data Export ({format.upper()}) ===")
    
    try:
        # Build query filters
        filters = QueryFilters(
            pfam_family=pfam_family,
            protein_id=protein_id,
            min_sequence_length=min_length,
            max_sequence_length=max_length,
            has_tim_barrel=tim_barrel if tim_barrel else None
        )
        
        # Query data
        if pfam_family:
            click.echo(f"Filtering by PFAM family: {pfam_family}")
            result = engine.search_by_pfam_family(pfam_family, True)
            proteins = result.isoforms
        elif protein_id:
            click.echo(f"Filtering by protein ID: {protein_id}")
            result = engine.get_protein_isoforms(protein_id)
            proteins = result.isoforms
        elif tim_barrel:
            click.echo("Filtering proteins with TIM barrel annotations")
            criteria = {'has_location': True}
            result = engine.search_by_tim_barrel_features(criteria)
            proteins = result.isoforms
        else:
            click.echo("Exporting all proteins with specified filters")
            result = engine.filter_proteins(filters, limit)
            proteins = result.isoforms
        
        if not proteins:
            click.echo("No proteins found matching criteria.")
            return
        
        click.echo(f"Found {len(proteins)} proteins to export")
        
        # Apply limit if specified
        if limit and len(proteins) > limit:
            proteins = proteins[:limit]
            click.echo(f"Limited to {limit} proteins")
        
        # Set up export options
        export_format = ExportFormat(format)
        options = ExportOptions(
            format=export_format,
            include_sequences=include_sequences,
            include_metadata=include_metadata,
            include_annotations=include_annotations,
            max_sequence_length=max_seq_length
        )
        
        # Export data
        exported_data = exporter.export_proteins(proteins, export_format, options)
        
        # Validate format
        if not exporter.validate_export_format(exported_data, export_format):
            click.echo("Warning: Exported data may not be valid for the specified format", err=True)
        
        # Output results
        if output:
            with open(output, 'w') as f:
                f.write(exported_data)
            click.echo(f"Data exported to: {output}")
        else:
            click.echo(exported_data)
        
        click.echo(f"\nExport completed: {len(proteins)} proteins in {format} format")
        
    except Exception as e:
        click.echo(f"Export failed: {e}", err=True)
        sys.exit(1)


@cli.command()
@click.pass_context
def status(ctx):
    """Show system status and configuration."""
    from .database.connection import get_database_manager
    from .query.engine import QueryEngine
    import requests
    
    config = ctx.obj['config']
    
    click.echo("=== Protein Data Collector Status ===")
    
    # Configuration info
    click.echo(f"\nConfiguration:")
    click.echo(f"  Database: {config.database.host}:{config.database.port}/{config.database.database}")
    click.echo(f"  InterPro API: {config.api.interpro_base_url}")
    click.echo(f"  UniProt API: {config.api.uniprot_base_url}")
    click.echo(f"  Max retries: {config.retry.max_retries}")
    click.echo(f"  Log level: {config.logging.level}")
    
    # Database connectivity check
    click.echo(f"\nDatabase Status:")
    try:
        db_manager = get_database_manager()
        if db_manager.test_connection():
            click.echo("  ✓ Database connection: OK")
            
            # Get data statistics
            engine = QueryEngine()
            stats = engine.get_summary_statistics()
            click.echo(f"  Data Summary:")
            click.echo(f"    PFAM families: {stats['pfam_families']}")
            click.echo(f"    Proteins: {stats['proteins']}")
            click.echo(f"    Isoforms: {stats['isoforms']}")
            click.echo(f"    TIM barrel annotations: {stats['tim_barrel_annotations']}")
            click.echo(f"    TIM barrel coverage: {stats['tim_barrel_coverage']:.1f}%")
            
            if stats['sequence_length']['avg']:
                click.echo(f"    Avg sequence length: {stats['sequence_length']['avg']:.0f} aa")
            if stats['exon_count']['avg']:
                click.echo(f"    Avg exon count: {stats['exon_count']['avg']:.1f}")
        else:
            click.echo("  ✗ Database connection: FAILED")
    except Exception as e:
        click.echo(f"  ✗ Database connection: ERROR - {e}")
    
    # API availability check
    click.echo(f"\nAPI Status:")
    
    # Check InterPro API
    try:
        response = requests.get(f"{config.api.interpro_base_url}entry/pfam/", 
                              timeout=5, params={'page_size': 1})
        if response.status_code == 200:
            click.echo("  ✓ InterPro API: OK")
        else:
            click.echo(f"  ✗ InterPro API: HTTP {response.status_code}")
    except Exception as e:
        click.echo(f"  ✗ InterPro API: ERROR - {e}")
    
    # Check UniProt API
    try:
        response = requests.get(f"{config.api.uniprot_base_url}uniprotkb/search", 
                              timeout=5, params={'query': 'organism_id:9606', 'size': 1})
        if response.status_code == 200:
            click.echo("  ✓ UniProt API: OK")
        else:
            click.echo(f"  ✗ UniProt API: HTTP {response.status_code}")
    except Exception as e:
        click.echo(f"  ✗ UniProt API: ERROR - {e}")

    else:
        click.echo("  UniProt: REST API only")


@cli.command()
@click.option('--progress-file', type=click.Path(exists=True), required=True,
              help='Progress file to check')
@click.pass_context
def progress(ctx, progress_file):
    """Check collection progress from progress file."""
    import json
    from datetime import datetime
    
    try:
        with open(progress_file, 'r') as f:
            progress_data = json.load(f)
        
        click.echo("=== Collection Progress ===")
        click.echo(f"Current phase: {progress_data.get('phase', 'unknown')}")
        click.echo(f"PFAM families collected: {progress_data.get('pfam_families_collected', 0)}")
        click.echo(f"InterPro proteins collected: {progress_data.get('interpro_proteins_collected', 0)}")
        click.echo(f"UniProt isoforms collected: {progress_data.get('uniprot_isoforms_collected', 0)}")
        click.echo(f"Total entities collected: {progress_data.get('total_entities_collected', 0)}")
        
        if progress_data.get('start_time'):
            start_time = datetime.fromisoformat(progress_data['start_time'])
            click.echo(f"Started: {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
        
        if progress_data.get('last_checkpoint'):
            checkpoint_time = datetime.fromisoformat(progress_data['last_checkpoint'])
            click.echo(f"Last checkpoint: {checkpoint_time.strftime('%Y-%m-%d %H:%M:%S')}")
        
        duration = progress_data.get('duration_seconds', 0)
        if duration > 0:
            hours = int(duration // 3600)
            minutes = int((duration % 3600) // 60)
            seconds = int(duration % 60)
            click.echo(f"Duration: {hours:02d}:{minutes:02d}:{seconds:02d}")
        
        if progress_data.get('errors'):
            click.echo(f"\nErrors encountered: {len(progress_data['errors'])}")
            for error in progress_data['errors'][-5:]:  # Show last 5 errors
                click.echo(f"  - {error}")
        
    except Exception as e:
        click.echo(f"Failed to read progress file: {e}", err=True)
        sys.exit(1)


@cli.command()
@click.option('--output', '-o', type=click.Path(), 
              help='Output configuration file path')
@click.pass_context
def config_export(ctx, output):
    """Export current configuration to file."""
    config = ctx.obj['config']
    
    if output:
        config.to_file(output)
        click.echo(f"Configuration exported to: {output}")
    else:
        # Print to stdout
        import json
        config_dict = {
            "database": {
                "host": config.database.host,
                "port": config.database.port,
                "database": config.database.database,
                "username": config.database.username
            },
            "api": {
                "interpro_base_url": config.api.interpro_base_url,
                "uniprot_base_url": config.api.uniprot_base_url
            },
            "retry": {
                "max_retries": config.retry.max_retries,
                "initial_delay": config.retry.initial_delay,
                "backoff_multiplier": config.retry.backoff_multiplier
            }
        }
        click.echo(json.dumps(config_dict, indent=2))


@cli.command()
@click.pass_context
def stats(ctx):
    """Show detailed database statistics."""
    from .query.engine import QueryEngine
    
    try:
        engine = QueryEngine()
        stats = engine.get_summary_statistics()
        
        click.echo("=== Database Statistics ===")
        click.echo(f"PFAM families: {stats['pfam_families']}")
        click.echo(f"Proteins: {stats['proteins']}")
        click.echo(f"Isoforms: {stats['isoforms']}")
        click.echo(f"TIM barrel annotations: {stats['tim_barrel_annotations']}")
        click.echo(f"TIM barrel coverage: {stats['tim_barrel_coverage']:.1f}%")
        
        click.echo(f"\nSequence Length Statistics:")
        if stats['sequence_length']['min']:
            click.echo(f"  Minimum: {stats['sequence_length']['min']} aa")
        if stats['sequence_length']['max']:
            click.echo(f"  Maximum: {stats['sequence_length']['max']} aa")
        if stats['sequence_length']['avg']:
            click.echo(f"  Average: {stats['sequence_length']['avg']:.1f} aa")
        
        click.echo(f"\nExon Count Statistics:")
        if stats['exon_count']['min']:
            click.echo(f"  Minimum: {stats['exon_count']['min']}")
        if stats['exon_count']['max']:
            click.echo(f"  Maximum: {stats['exon_count']['max']}")
        if stats['exon_count']['avg']:
            click.echo(f"  Average: {stats['exon_count']['avg']:.1f}")
        
    except Exception as e:
        click.echo(f"Failed to get statistics: {e}", err=True)
        sys.exit(1)


@cli.command()
@click.option('--host', default='0.0.0.0', help='Host to bind to')
@click.option('--port', default=8000, type=int, help='Port to bind to')
@click.option('--reload', is_flag=True, help='Enable auto-reload for development')
@click.option('--workers', default=1, type=int, help='Number of worker processes')
@click.pass_context
def serve(ctx, host, port, reload, workers):
    """Start the REST API server."""
    from .server import run_server
    
    click.echo(f"=== Starting Protein Data Collector API Server ===")
    click.echo(f"Host: {host}")
    click.echo(f"Port: {port}")
    click.echo(f"Reload: {reload}")
    click.echo(f"Workers: {workers}")
    click.echo(f"API Documentation: http://{host}:{port}/docs")
    click.echo(f"Health Check: http://{host}:{port}/health")
    
    try:
        run_server(host=host, port=port, reload=reload, workers=workers)
    except KeyboardInterrupt:
        click.echo("\nServer stopped by user.")
    except Exception as e:
        click.echo(f"Server failed to start: {e}", err=True)
        sys.exit(1)


def main():
    """Main entry point for the CLI."""
    try:
        cli()
    except KeyboardInterrupt:
        click.echo("\nOperation cancelled by user.", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


if __name__ == '__main__':
    main()