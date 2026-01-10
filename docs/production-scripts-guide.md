# Production Scripts Guide

## Overview
This document describes the production-ready scripts for the TIM Barrel Protein Data Collector project.

## Production Scripts

### 1. `collect_tim_barrel_entries.py`
**Purpose**: Collects TIM barrel domain entries from InterPro and Pfam databases
**Usage**: Initial setup to populate `tim_barrel_entries` table
**Status**: Production-ready, well-documented

### 2. `collect_human_proteins.py`
**Purpose**: Collects human proteins associated with TIM barrel domains
**Usage**: Populates `interpro_proteins` table with protein-domain associations
**Status**: Production-ready, well-documented

### 3. `process_all_407_proteins.py`
**Purpose**: Main production collector for complete isoform data
**Features**:
- Processes all 407 proteins with complete isoform data
- Collects sequences, domain boundaries, and Ensembl references
- 100% success rate achieved
- Real data only - no hardcoded fallbacks
**Usage**: Final production run for complete database population
**Status**: Production-ready, successfully tested

### 4. `direct_uniprot_api_collector.py`
**Purpose**: Direct UniProt API collector for individual protein processing
**Features**:
- Direct UniProt REST API integration
- Complete isoform data extraction
- Domain boundary detection
- Ensembl reference mapping
**Usage**: Individual protein processing or debugging
**Status**: Production-ready, core component

### 5. `create_sqlite_simple.py`
**Purpose**: Database schema creation and initialization
**Usage**: Sets up the SQLite database with proper schema
**Status**: Production-ready

## Database Schema
- **tim_barrel_entries**: TIM barrel domain definitions
- **interpro_proteins**: Protein-domain associations
- **proteins**: Complete protein data with isoform information

## Data Quality
- **407 proteins** total in database
- **494 isoforms** with complete sequences
- **100% TIM barrel boundary coverage**
- **38 multi-isoform proteins** for alternative splicing research

## Usage Workflow
1. Run `create_sqlite_simple.py` to create database
2. Run `collect_tim_barrel_entries.py` to populate domain entries
3. Run `collect_human_proteins.py` to populate protein associations
4. Run `process_all_407_proteins.py` for complete data collection

## Development Scripts
Development and testing scripts are located in `scripts/temp/` directory.