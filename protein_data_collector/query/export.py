"""
Data export functionality for standard bioinformatics formats.

This module provides export methods for FASTA, JSON, and CSV formats with
format validation and streaming support for large datasets.
"""

import json
import csv
import logging
from typing import List, Dict, Any, Iterator, Optional, Union, TextIO
from io import StringIO
from dataclasses import dataclass
from enum import Enum

logger = logging.getLogger(__name__)


class ExportFormat(Enum):
    """Supported export formats."""
    FASTA = "fasta"
    JSON = "json"
    CSV = "csv"


@dataclass
class ExportOptions:
    """Configuration options for data export."""
    format: ExportFormat
    include_metadata: bool = True
    include_sequences: bool = True
    include_annotations: bool = True
    max_sequence_length: Optional[int] = None
    pretty_json: bool = True
    csv_delimiter: str = ","
    fasta_line_length: int = 80


class DataExporter:
    """
    Data exporter for protein data in standard bioinformatics formats.
    
    Supports FASTA, JSON, and CSV formats with format validation and
    streaming capabilities for large datasets.
    """
    
    def __init__(self):
        """Initialize the data exporter."""
        self.logger = logging.getLogger(__name__)
    
    def export_proteins(
        self, 
        proteins: List[Dict[str, Any]], 
        format: ExportFormat, 
        options: Optional[ExportOptions] = None
    ) -> str:
        """
        Export protein data in the specified format.
        
        Args:
            proteins: List of protein dictionaries to export
            format: Export format (FASTA, JSON, or CSV)
            options: Export configuration options
            
        Returns:
            Formatted string containing exported data
            
        Requirements: 5.6
        """
        if options is None:
            options = ExportOptions(format=format)
        
        self.logger.info(f"Exporting {len(proteins)} proteins in {format.value} format")
        
        if format == ExportFormat.FASTA:
            return self._export_fasta(proteins, options)
        elif format == ExportFormat.JSON:
            return self._export_json(proteins, options)
        elif format == ExportFormat.CSV:
            return self._export_csv(proteins, options)
        else:
            raise ValueError(f"Unsupported export format: {format}")
    
    def export_proteins_stream(
        self, 
        proteins: Iterator[Dict[str, Any]], 
        output: TextIO,
        format: ExportFormat, 
        options: Optional[ExportOptions] = None
    ) -> int:
        """
        Export protein data as a stream for large datasets.
        
        Args:
            proteins: Iterator of protein dictionaries to export
            output: Text stream to write to
            format: Export format (FASTA, JSON, or CSV)
            options: Export configuration options
            
        Returns:
            Number of proteins exported
            
        Requirements: 5.6
        """
        if options is None:
            options = ExportOptions(format=format)
        
        self.logger.info(f"Starting streaming export in {format.value} format")
        
        if format == ExportFormat.FASTA:
            return self._export_fasta_stream(proteins, output, options)
        elif format == ExportFormat.JSON:
            return self._export_json_stream(proteins, output, options)
        elif format == ExportFormat.CSV:
            return self._export_csv_stream(proteins, output, options)
        else:
            raise ValueError(f"Unsupported export format: {format}")
    
    def validate_export_format(self, data: str, format: ExportFormat) -> bool:
        """
        Validate that exported data meets format specifications.
        
        Args:
            data: Exported data string to validate
            format: Expected format
            
        Returns:
            True if data is valid for the format
            
        Requirements: 5.6
        """
        try:
            if format == ExportFormat.FASTA:
                return self._validate_fasta_format(data)
            elif format == ExportFormat.JSON:
                return self._validate_json_format(data)
            elif format == ExportFormat.CSV:
                return self._validate_csv_format(data)
            else:
                return False
        except Exception as e:
            self.logger.error(f"Format validation error: {e}")
            return False
    
    def _export_fasta(self, proteins: List[Dict[str, Any]], options: ExportOptions) -> str:
        """Export proteins in FASTA format."""
        output = StringIO()
        
        for protein in proteins:
            if not options.include_sequences or not protein.get('sequence'):
                continue
            
            # Create FASTA header
            header_parts = []
            
            # Primary identifier
            protein_id = protein.get('isoform_id') or protein.get('uniprot_id', 'unknown')
            header_parts.append(protein_id)
            
            # Add description/name
            if protein.get('name'):
                header_parts.append(protein['name'])
            elif protein.get('description'):
                header_parts.append(protein['description'])
            
            # Add organism if available
            if options.include_metadata and protein.get('organism'):
                header_parts.append(f"OS={protein['organism']}")
            
            # Add sequence length
            if options.include_metadata and protein.get('sequence_length'):
                header_parts.append(f"Length={protein['sequence_length']}")
            
            # Add TIM barrel annotation if available
            if options.include_annotations and protein.get('tim_barrel_location'):
                tim_loc = protein['tim_barrel_location']
                if isinstance(tim_loc, dict) and 'start' in tim_loc and 'end' in tim_loc:
                    header_parts.append(f"TIM_barrel={tim_loc['start']}-{tim_loc['end']}")
            
            # Write header
            header = " | ".join(header_parts)
            output.write(f">{header}\n")
            
            # Write sequence with line wrapping
            sequence = protein['sequence']
            if options.max_sequence_length:
                sequence = sequence[:options.max_sequence_length]
            
            line_length = options.fasta_line_length
            for i in range(0, len(sequence), line_length):
                output.write(sequence[i:i + line_length] + "\n")
        
        return output.getvalue()
    
    def _export_json(self, proteins: List[Dict[str, Any]], options: ExportOptions) -> str:
        """Export proteins in JSON format."""
        export_data = []
        
        for protein in proteins:
            protein_data = {}
            
            # Always include basic identifiers
            for key in ['isoform_id', 'uniprot_id', 'parent_protein_id']:
                if key in protein:
                    protein_data[key] = protein[key]
            
            # Include sequences if requested
            if options.include_sequences:
                for key in ['sequence', 'sequence_length']:
                    if key in protein:
                        value = protein[key]
                        if key == 'sequence' and options.max_sequence_length:
                            value = value[:options.max_sequence_length]
                        protein_data[key] = value
            
            # Include metadata if requested
            if options.include_metadata:
                for key in ['name', 'description', 'organism', 'created_at']:
                    if key in protein:
                        protein_data[key] = protein[key]
            
            # Include annotations if requested
            if options.include_annotations:
                for key in ['exon_annotations', 'exon_count', 'tim_barrel_location', 'pfam_accession']:
                    if key in protein:
                        protein_data[key] = protein[key]
            
            export_data.append(protein_data)
        
        # Export as JSON
        if options.pretty_json:
            return json.dumps(export_data, indent=2, ensure_ascii=False)
        else:
            return json.dumps(export_data, ensure_ascii=False)
    
    def _export_csv(self, proteins: List[Dict[str, Any]], options: ExportOptions) -> str:
        """Export proteins in CSV format."""
        if not proteins:
            return ""
        
        output = StringIO()
        
        # Determine columns based on options and available data
        columns = ['isoform_id', 'uniprot_id', 'parent_protein_id']
        
        if options.include_sequences:
            columns.extend(['sequence_length'])
            if any(p.get('sequence') for p in proteins):
                columns.append('sequence')
        
        if options.include_metadata:
            columns.extend(['name', 'description', 'organism', 'created_at'])
        
        if options.include_annotations:
            columns.extend(['pfam_accession', 'exon_count'])
            if any(p.get('tim_barrel_location') for p in proteins):
                columns.extend(['tim_barrel_start', 'tim_barrel_end', 'tim_barrel_confidence'])
        
        # Remove columns that don't exist in any protein
        available_columns = []
        for col in columns:
            if col.startswith('tim_barrel_'):
                # Special handling for TIM barrel location fields
                if any(p.get('tim_barrel_location') for p in proteins):
                    available_columns.append(col)
            elif any(col in p for p in proteins):
                available_columns.append(col)
        
        writer = csv.DictWriter(
            output, 
            fieldnames=available_columns, 
            delimiter=options.csv_delimiter,
            extrasaction='ignore'
        )
        writer.writeheader()
        
        for protein in proteins:
            row = {}
            
            # Copy basic fields
            for col in available_columns:
                if col.startswith('tim_barrel_'):
                    # Extract TIM barrel location fields
                    tim_loc = protein.get('tim_barrel_location')
                    if isinstance(tim_loc, dict):
                        if col == 'tim_barrel_start':
                            row[col] = tim_loc.get('start')
                        elif col == 'tim_barrel_end':
                            row[col] = tim_loc.get('end')
                        elif col == 'tim_barrel_confidence':
                            row[col] = tim_loc.get('confidence')
                elif col == 'sequence' and options.max_sequence_length:
                    sequence = protein.get('sequence', '')
                    row[col] = sequence[:options.max_sequence_length] if sequence else ''
                else:
                    row[col] = protein.get(col, '')
            
            writer.writerow(row)
        
        return output.getvalue()
    
    def _export_fasta_stream(
        self, 
        proteins: Iterator[Dict[str, Any]], 
        output: TextIO, 
        options: ExportOptions
    ) -> int:
        """Stream export proteins in FASTA format."""
        count = 0
        
        for protein in proteins:
            if not options.include_sequences or not protein.get('sequence'):
                continue
            
            # Create FASTA header
            header_parts = []
            
            # Primary identifier
            protein_id = protein.get('isoform_id') or protein.get('uniprot_id', 'unknown')
            header_parts.append(protein_id)
            
            # Add description/name
            if protein.get('name'):
                header_parts.append(protein['name'])
            elif protein.get('description'):
                header_parts.append(protein['description'])
            
            # Add organism if available
            if options.include_metadata and protein.get('organism'):
                header_parts.append(f"OS={protein['organism']}")
            
            # Add sequence length
            if options.include_metadata and protein.get('sequence_length'):
                header_parts.append(f"Length={protein['sequence_length']}")
            
            # Add TIM barrel annotation if available
            if options.include_annotations and protein.get('tim_barrel_location'):
                tim_loc = protein['tim_barrel_location']
                if isinstance(tim_loc, dict) and 'start' in tim_loc and 'end' in tim_loc:
                    header_parts.append(f"TIM_barrel={tim_loc['start']}-{tim_loc['end']}")
            
            # Write header
            header = " | ".join(header_parts)
            output.write(f">{header}\n")
            
            # Write sequence with line wrapping
            sequence = protein['sequence']
            if options.max_sequence_length:
                sequence = sequence[:options.max_sequence_length]
            
            line_length = options.fasta_line_length
            for i in range(0, len(sequence), line_length):
                output.write(sequence[i:i + line_length] + "\n")
            
            count += 1
        
        return count
    
    def _export_json_stream(
        self, 
        proteins: Iterator[Dict[str, Any]], 
        output: TextIO, 
        options: ExportOptions
    ) -> int:
        """Stream export proteins in JSON format."""
        count = 0
        output.write("[\n")
        
        first = True
        for protein in proteins:
            if not first:
                output.write(",\n")
            first = False
            
            protein_data = {}
            
            # Always include basic identifiers
            for key in ['isoform_id', 'uniprot_id', 'parent_protein_id']:
                if key in protein:
                    protein_data[key] = protein[key]
            
            # Include sequences if requested
            if options.include_sequences:
                for key in ['sequence', 'sequence_length']:
                    if key in protein:
                        value = protein[key]
                        if key == 'sequence' and options.max_sequence_length:
                            value = value[:options.max_sequence_length]
                        protein_data[key] = value
            
            # Include metadata if requested
            if options.include_metadata:
                for key in ['name', 'description', 'organism', 'created_at']:
                    if key in protein:
                        protein_data[key] = protein[key]
            
            # Include annotations if requested
            if options.include_annotations:
                for key in ['exon_annotations', 'exon_count', 'tim_barrel_location', 'pfam_accession']:
                    if key in protein:
                        protein_data[key] = protein[key]
            
            # Write JSON object
            if options.pretty_json:
                json_str = json.dumps(protein_data, indent=2, ensure_ascii=False)
                # Indent the entire object
                indented_lines = ["  " + line for line in json_str.split("\n")]
                output.write("\n".join(indented_lines))
            else:
                output.write(json.dumps(protein_data, ensure_ascii=False))
            
            count += 1
        
        output.write("\n]\n")
        return count
    
    def _export_csv_stream(
        self, 
        proteins: Iterator[Dict[str, Any]], 
        output: TextIO, 
        options: ExportOptions
    ) -> int:
        """Stream export proteins in CSV format."""
        count = 0
        writer = None
        
        for protein in proteins:
            if writer is None:
                # Initialize writer with first protein to determine columns
                columns = ['isoform_id', 'uniprot_id', 'parent_protein_id']
                
                if options.include_sequences:
                    columns.extend(['sequence_length'])
                    if protein.get('sequence'):
                        columns.append('sequence')
                
                if options.include_metadata:
                    columns.extend(['name', 'description', 'organism', 'created_at'])
                
                if options.include_annotations:
                    columns.extend(['pfam_accession', 'exon_count'])
                    if protein.get('tim_barrel_location'):
                        columns.extend(['tim_barrel_start', 'tim_barrel_end', 'tim_barrel_confidence'])
                
                # Filter to available columns
                available_columns = []
                for col in columns:
                    if col.startswith('tim_barrel_'):
                        if protein.get('tim_barrel_location'):
                            available_columns.append(col)
                    elif col in protein:
                        available_columns.append(col)
                
                writer = csv.DictWriter(
                    output, 
                    fieldnames=available_columns, 
                    delimiter=options.csv_delimiter,
                    extrasaction='ignore'
                )
                writer.writeheader()
            
            # Write row
            row = {}
            for col in writer.fieldnames:
                if col.startswith('tim_barrel_'):
                    # Extract TIM barrel location fields
                    tim_loc = protein.get('tim_barrel_location')
                    if isinstance(tim_loc, dict):
                        if col == 'tim_barrel_start':
                            row[col] = tim_loc.get('start')
                        elif col == 'tim_barrel_end':
                            row[col] = tim_loc.get('end')
                        elif col == 'tim_barrel_confidence':
                            row[col] = tim_loc.get('confidence')
                elif col == 'sequence' and options.max_sequence_length:
                    sequence = protein.get('sequence', '')
                    row[col] = sequence[:options.max_sequence_length] if sequence else ''
                else:
                    row[col] = protein.get(col, '')
            
            writer.writerow(row)
            count += 1
        
        return count
    
    def _validate_fasta_format(self, data: str) -> bool:
        """Validate FASTA format."""
        if not data.strip():
            return True  # Empty data is valid
        
        lines = data.strip().split('\n')
        expecting_header = True
        
        for line in lines:
            line = line.strip()
            if not line:
                continue
            
            if line.startswith('>'):
                if not expecting_header:
                    # Found header when expecting sequence - this is okay
                    pass
                expecting_header = False
            else:
                if expecting_header:
                    return False  # Found sequence without header
                # Validate sequence contains only valid characters
                if not all(c.upper() in 'ACDEFGHIKLMNPQRSTVWYX-' for c in line):
                    return False
        
        return True
    
    def _validate_json_format(self, data: str) -> bool:
        """Validate JSON format."""
        try:
            parsed = json.loads(data)
            # Should be a list of objects
            if not isinstance(parsed, list):
                return False
            
            for item in parsed:
                if not isinstance(item, dict):
                    return False
            
            return True
        except (json.JSONDecodeError, TypeError):
            return False
    
    def _validate_csv_format(self, data: str) -> bool:
        """Validate CSV format."""
        if not data.strip():
            return True  # Empty data is valid
        
        try:
            # Try to parse as CSV
            reader = csv.reader(StringIO(data))
            rows = list(reader)
            
            if not rows:
                return True  # Empty CSV is valid
            
            # Check that all rows have the same number of columns
            header_length = len(rows[0])
            for row in rows[1:]:  # Skip header
                if len(row) != header_length:
                    return False
            
            return True
        except (csv.Error, UnicodeDecodeError):
            return False