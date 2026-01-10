# UniProt Fields Summary - Complete Field Inventory

## Total Field Count: **250+ Fields**

UniProt provides extensive protein annotation data through their REST API. Here's a comprehensive breakdown of all available fields organized by category:

## 📊 Field Categories and Counts

### 1. **Names & Taxonomy** (14 fields)
- `accession` - Entry accession
- `id` - Entry Name  
- `gene_names` - Gene Names (all)
- `gene_primary` - Gene Names (primary)
- `gene_synonym` - Gene Names (synonym)
- `gene_oln` - Gene Names (ordered locus)
- `gene_orf` - Gene Names (ORF)
- `organism_name` - Organism
- `organism_id` - Organism ID
- `protein_name` - Protein names
- `xref_proteomes` - Proteomes
- `lineage` - Taxonomic lineage
- `lineage_ids` - Taxonomic lineage (IDs)
- `virus_hosts` - Virus hosts

### 2. **Sequences** (19 fields)
- `cc_alternative_products` - Alternative products
- `ft_var_seq` - Alternative sequence
- `error_gmodel_pred` - Erroneous gene model prediction
- `fragment` - Fragment
- `organelle` - Gene encoded by
- `length` - Length
- `mass` - Mass
- `cc_mass_spectrometry` - Mass spectrometry
- `ft_variant` - Natural variant
- `ft_non_cons` - Non-adjacent residues
- `ft_non_std` - Non-standard residue
- `ft_non_ter` - Non-terminal residue
- `cc_polymorphism` - Polymorphism
- `cc_rna_editing` - RNA editing
- `sequence` - Sequence
- `cc_sequence_caution` - Sequence caution
- `ft_conflict` - Sequence conflict
- `ft_unsure` - Sequence uncertainty
- `sequence_version` - Sequence version

### 3. **Function** (15 fields)
- `absorption` - Absorption
- `ft_act_site` - Active site
- `cc_activity_regulation` - Activity regulation
- `ft_binding` - Binding site
- `cc_catalytic_activity` - Catalytic activity
- `cc_cofactor` - Cofactor
- `ft_dna_bind` - DNA binding
- `ec` - EC number
- `cc_function` - Function [CC]
- `kinetics` - Kinetics
- `cc_pathway` - Pathway
- `ph_dependence` - pH dependence
- `redox_potential` - Redox potential
- `rhea` - Rhea ID
- `ft_site` - Site
- `temp_dependence` - Temperature dependence

### 4. **Miscellaneous** (11 fields)
- `annotation_score` - Annotation
- `cc_caution` - Caution
- `comment_count` - Comment Count
- `feature_count` - Features
- `keywordid` - Keyword ID
- `keyword` - Keywords
- `cc_miscellaneous` - Miscellaneous [CC]
- `protein_existence` - Protein existence
- `reviewed` - Reviewed
- `tools` - Tools
- `uniparc_id` - UniParc

### 5. **Interaction** (2 fields)
- `cc_interaction` - Interacts with
- `cc_subunit` - Subunit structure[CC]

### 6. **Expression** (3 fields)
- `cc_developmental_stage` - Developmental stage
- `cc_induction` - Induction
- `cc_tissue_specificity` - Tissue specificity

### 7. **Gene Ontology (GO)** (5 fields)
- `go_p` - Gene ontology (biological process)
- `go_c` - Gene ontology (cellular component)
- `go` - Gene ontology (GO)
- `go_f` - Gene ontology (molecular function)
- `go_id` - Gene ontology IDs

### 8. **Pathology & Biotech** (7 fields)
- `cc_allergen` - Allergenic properties
- `cc_biotechnology` - Biotechnological use
- `cc_disruption_phenotype` - Disruption phenotype
- `cc_disease` - Involvement in disease
- `ft_mutagen` - Mutagenesis
- `cc_pharmaceutical` - Pharmaceutical use
- `cc_toxic_dose` - Toxic dose

### 9. **Subcellular Location** (4 fields)
- `ft_intramem` - Intramembrane
- `cc_subcellular_location` - Subcellular location[CC]
- `ft_topo_dom` - Topological domain
- `ft_transmem` - Transmembrane

### 10. **PTM / Processing** (12 fields)
- `ft_chain` - Chain
- `ft_crosslnk` - Cross-link
- `ft_disulfid` - Disulfide bond
- `ft_carbohyd` - Glycosylation
- `ft_init_met` - Initiator methionine
- `ft_lipid` - Lipidation
- `ft_mod_res` - Modified residue
- `ft_peptide` - Peptide
- `cc_ptm` - Post-translational modification
- `ft_propep` - Propeptide
- `ft_signal` - Signal peptide
- `ft_transit` - Transit peptide

### 11. **Structure** (4 fields)
- `structure_3d` - 3D
- `ft_strand` - Beta strand
- `ft_helix` - Helix
- `ft_turn` - Turn

### 12. **Publications** (1 field)
- `lit_pubmed_id` - PubMed ID

### 13. **Date Information** (4 fields)
- `date_created` - Date of creation
- `date_modified` - Date of last modification
- `date_sequence_modified` - Date of last sequence modification
- `version` - Entry version

### 14. **Family & Domains** (9 fields)
- `ft_coiled` - Coiled coil
- `ft_compbias` - Compositional bias
- `cc_domain` - Domain[CC]
- `ft_domain` - Domain[FT]
- `ft_motif` - Motif
- `protein_families` - Protein families
- `ft_region` - Region
- `ft_repeat` - Repeat
- `ft_zn_fing` - Zinc finger

## 🔗 Cross-Reference Databases (150+ fields)

### **Sequence Databases** (4 fields)
- `xref_ccds` - CCDS
- `xref_embl` - EMBL
- `xref_pir` - PIR
- `xref_refseq` - RefSeq

### **3D Structure Databases** (7 fields)
- `xref_alphafolddb` - AlphaFoldDB
- `xref_bmrb` - BMRB
- `xref_pcddb` - PCDDB
- `xref_pdb` - PDB
- `xref_pdbsum` - PDBsum
- `xref_sasbdb` - SASBDB
- `xref_smr` - SMR

### **Protein-Protein Interaction Databases** (8 fields)
- `xref_biogrid` - BioGRID
- `xref_corum` - CORUM
- `xref_complexportal` - ComplexPortal
- `xref_dip` - DIP
- `xref_elm` - ELM
- `xref_intact` - IntAct
- `xref_mint` - MINT
- `xref_string` - STRING

### **Chemistry Databases** (6 fields)
- `xref_bindingdb` - BindingDB
- `xref_chembl` - ChEMBL
- `xref_drugbank` - DrugBank
- `xref_drugcentral` - DrugCentral
- `xref_guidetopharmacology` - GuidetoPHARMACOLOGY
- `xref_swisslipids` - SwissLipids

### **Protein Family/Group Databases** (12 fields)
- `xref_allergome` - Allergome
- `xref_cazy` - CAZy
- `xref_clae` - CLAE
- `xref_esther` - ESTHER
- `xref_imgt_gene-db` - IMGT_GENE-DB
- `xref_merops` - MEROPS
- `xref_moondb` - MoonDB
- `xref_moonprot` - MoonProt
- `xref_peroxibase` - PeroxiBase
- `xref_rebase` - REBASE
- `xref_tcdb` - TCDB
- `xref_unilectin` - UniLectin

### **PTM Databases** (9 fields)
- `xref_carbonyldb` - CarbonylDB
- `xref_depod` - DEPOD
- `xref_glycosmos` - GlyCosmos
- `xref_glyconnect` - GlyConnect
- `xref_glygen` - GlyGen
- `xref_metosite` - MetOSite
- `xref_phosphositeplus` - PhosphoSitePlus
- `xref_swisspalm` - SwissPalm
- `xref_iptmnet` - iPTMnet

### **Genetic Variation Databases** (3 fields)
- `xref_biomuta` - BioMuta
- `xref_dmdm` - DMDM
- `xref_dbsnp` - dbSNP

### **2D Gel Databases** (7 fields)
- `xref_compluyeast-2dpage` - COMPLUYEAST-2DPAGE
- `xref_dosac-cobs-2dpage` - DOSAC-COBS-2DPAGE
- `xref_ogp` - OGP
- `xref_reproduction-2dpage` - REPRODUCTION-2DPAGE
- `xref_swiss-2dpage` - SWISS-2DPAGE
- `xref_ucd-2dpage` - UCD-2DPAGE
- `xref_world-2dpage` - World-2DPAGE

### **Proteomic Databases** (11 fields)
- `xref_cptac` - CPTAC
- `xref_epd` - EPD
- `xref_massive` - MassIVE
- `xref_maxqb` - MaxQB
- `xref_pride` - PRIDE
- `xref_paxdb` - PaxDb
- `xref_peptideatlas` - PeptideAtlas
- `xref_promex` - ProMEX
- `xref_proteomicsdb` - ProteomicsDB
- `xref_topdownproteomics` - TopDownProteomics
- `xref_jpost` - jPOST

### **Protocols And Materials Databases** (4 fields)
- `xref_abcd` - ABCD
- `xref_antibodypedia` - Antibodypedia
- `xref_cptc` - CPTC
- `xref_dnasu` - DNASU

### **Genome Annotation Databases** (15 fields)
- `xref_ensembl` - Ensembl
- `xref_ensemblbacteria` - EnsemblBacteria
- `xref_ensemblfungi` - EnsemblFungi
- `xref_ensemblmetazoa` - EnsemblMetazoa
- `xref_ensemblplants` - EnsemblPlants
- `xref_ensemblprotists` - EnsemblProtists
- `xref_geneid` - GeneID
- `xref_gramene` - Gramene
- `xref_kegg` - KEGG
- `xref_mane-select` - MANE-Select
- `xref_patric` - PATRIC
- `xref_ucsc` - UCSC
- `xref_vectorbase` - VectorBase
- `xref_wbparasite` - WBParaSite
- `xref_wbparasitetranscriptprotein` - WBParaSiteTranscriptProtein

### **Organism-Specific Databases** (36 fields)
- `xref_agr` - AGR
- `xref_arachnoserver` - ArachnoServer
- `xref_araport` - Araport
- `xref_cgd` - CGD
- `xref_ctd` - CTD
- `xref_conoserver` - ConoServer
- `xref_disgenet` - DisGeNET
- `xref_echobase` - EchoBASE
- `xref_flybase` - FlyBase
- `xref_genecards` - GeneCards
- `xref_genereviews` - GeneReviews
- `xref_hgnc` - HGNC
- `xref_hpa` - HPA
- `xref_legiolist` - LegioList
- `xref_leproma` - Leproma
- `xref_mgi` - MGI
- `xref_mim` - MIM
- `xref_maizegdb` - MaizeGDB
- `xref_malacards` - MalaCards
- `xref_niagads` - NIAGADS
- `xref_opentargets` - OpenTargets
- `xref_orphanet` - Orphanet
- `xref_pharmgkb` - PharmGKB
- `xref_pombase` - PomBase
- `xref_pseudocap` - PseudoCAP
- `xref_rgd` - RGD
- `xref_sgd` - SGD
- `xref_tair` - TAIR
- `xref_tuberculist` - TubercuList
- `xref_veupathdb` - VEuPathDB
- `xref_vgnc` - VGNC
- `xref_wormbase` - WormBase
- `xref_xenbase` - Xenbase
- `xref_zfin` - ZFIN
- `xref_dictybase` - dictyBase
- `xref_euhcvdb` - euHCVdb
- `xref_nextprot` - neXtProt

### **Phylogenomic Databases** (9 fields)
- `xref_genetree` - GeneTree
- `xref_hogenom` - HOGENOM
- `xref_inparanoid` - InParanoid
- `xref_ko` - KO
- `xref_oma` - OMA
- `xref_orthodb` - OrthoDB
- `xref_phylomedb` - PhylomeDB
- `xref_treefam` - TreeFam
- `xref_eggnog` - eggNOG

### **Enzyme And Pathway Databases** (9 fields)
- `xref_brenda` - BRENDA
- `xref_biocyc` - BioCyc
- `xref_pathwaycommons` - PathwayCommons
- `xref_plantreactome` - PlantReactome
- `xref_reactome` - Reactome
- `xref_sabio-rk` - SABIO-RK
- `xref_signor` - SIGNOR
- `xref_signalink` - SignaLink
- `xref_unipathway` - UniPathway

### **Miscellaneous Databases** (9 fields)
- `xref_biogrid-orcs` - BioGRID-ORCS
- `xref_chitars` - ChiTaRS
- `xref_evolutionarytrace` - EvolutionaryTrace
- `xref_genewiki` - GeneWiki
- `xref_genomernai` - GenomeRNAi
- `xref_phi-base` - PHI-base
- `xref_pro` - PRO
- `xref_pharos` - Pharos
- `xref_rnact` - RNAct

### **Gene Expression Databases** (5 fields)
- `xref_bgee` - Bgee
- `xref_cleanex` - CleanEx
- `xref_collectf` - CollecTF
- `xref_expressionatlas` - ExpressionAtlas
- `xref_genevisible` - Genevisible

### **Family And Domain Databases** (16 fields)
- `xref_cdd` - CDD
- `xref_disprot` - DisProt
- `xref_gene3d` - Gene3D
- `xref_hamap` - HAMAP
- `xref_ideal` - IDEAL
- `xref_interpro` - InterPro
- `xref_panther` - PANTHER
- `xref_pirsf` - PIRSF
- `xref_prints` - PRINTS
- `xref_prosite` - PROSITE
- `xref_pfam` - Pfam
- `xref_prodom` - ProDom
- `xref_sfld` - SFLD
- `xref_smart` - SMART
- `xref_supfam` - SUPFAM
- `xref_tigrfams` - TIGRFAMs

## 📈 Field Population Statistics

**Note**: Not all proteins will have data for all fields. Field population varies by:

- **Reviewed vs Unreviewed**: Swiss-Prot (reviewed) entries have more complete annotations
- **Organism**: Human proteins typically have more comprehensive annotations
- **Research Interest**: Well-studied proteins have more experimental data
- **Database Coverage**: Some cross-references are organism or domain-specific

### Typical Field Population Rates:
- **Core Fields** (sequence, name, organism): ~100%
- **Functional Annotations**: 60-90% for reviewed entries
- **Structural Features**: 30-70% depending on experimental data
- **Cross-References**: Highly variable (5-95% depending on database)
- **Experimental Data**: 10-50% for specialized databases

## 🎯 Recommendations for TIM Barrel Research

### **Essential Fields** (Always collect):
1. Core identifiers and sequence data
2. Structural features (helix, strand, turn, domain)
3. Alternative products/isoforms
4. Gene Ontology terms
5. PDB and AlphaFold references

### **High-Value Fields** (Collect when available):
1. Functional annotations (catalytic activity, pathway)
2. InterPro and Pfam domain annotations
3. Tissue specificity and expression data
4. Disease associations
5. Protein-protein interactions

### **Research-Specific Fields** (Collect for specialized analysis):
1. Mass spectrometry data
2. Post-translational modifications
3. Genetic variations
4. Phylogenomic data
5. Expression databases

This comprehensive schema allows you to collect the complete UniProt dataset for your TIM barrel proteins, enabling sophisticated bioinformatics analyses and research insights.