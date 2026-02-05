import logging
import os
import gzip
import requests
import pandas as pd
from io import BytesIO
from gzip import BadGzipFile

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

# Constants
GDC_BASE = "https://api.gdc.cancer.gov"
PROJECT_ID = "TCGA-GBM"

# H1_GENES = [
#     "H1-1",
#     "H1-2", 
#     "H1-3",
#     "H1-5",
#     "H1-10",
#     "H1-4"
# ]

H1_GENES = [
   "TP53",      # Very common in GBM
    "PTEN",      # Common in GBM  
    "EGFR",      # Frequently amplified in GBM
    "IDH1"
    ]


OUTPUT_DIR = "data/processed/gbm"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# -------------------------------------------------
# Helper: GDC file search
# -------------------------------------------------
def search_gdc_files(data_type, data_format=None):
    filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": [PROJECT_ID]}},
            {"op": "in", "content": {"field": "data_type", "value": [data_type]}}
        ]
    }
    
    if data_format:
        filters["content"].append(
            {"op": "in", "content": {"field": "data_format", "value": [data_format]}}
        )

    params = {
        "filters": filters,
        "fields": "file_id,file_name,data_format,data_category",
        "format": "JSON",
        "size": 1000
    }

    try:
        r = requests.post(f"{GDC_BASE}/files", json=params, timeout=30)
        r.raise_for_status()
        return r.json()["data"]["hits"]
    except requests.exceptions.RequestException as e:
        logging.error(f"Error searching GDC files: {e}")
        return []

# -------------------------------------------------
# Download file by UUID
# -------------------------------------------------
def download_file(file_id):
    url = f"{GDC_BASE}/data/{file_id}"
    try:
        r = requests.get(url, stream=True, timeout=60)
        r.raise_for_status()
        return r.content
    except requests.exceptions.RequestException as e:
        logging.error(f"Error downloading file {file_id}: {e}")
        return None

# -------------------------------------------------
# Get gene symbol mapping
# -------------------------------------------------
def get_gene_symbol_mapping():
    """Get mapping from ENSEMBL IDs to gene symbols"""
    # Try to download a mapping file or use a local resource
    mapping_urls = [
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz",
        "https://www.gencodegenes.org/human/release_44/gencode.v44.annotation.gtf.gz"
    ]
    
    # Common H1 gene ENSEMBL IDs (these might vary)
    # h1_ensembl_map = {
    #     # These are common mappings - might need adjustment
    #     "H1-1": ["ENSG00000188153", "ENSG00000275373"],
    #     "H1-2": ["ENSG00000168496", "ENSG00000275600"],
    #     "H1-3": ["ENSG00000168484", "ENSG00000275023"],
    #     "H1-4": ["ENSG00000168495", "ENSG00000274286"],
    #     "H1-5": ["ENSG00000168490", "ENSG00000273748"],
    #     "H1-10": ["ENSG00000273125", "ENSG00000274444"]
    # }
    h1_ensembl_map = {
        # These are common mappings - might need adjustment
    "TP53":["ENSG00000141510"] ,      # Very common in GBM
    "PTEN":["ENSG00000171862"] ,      # Common in GBM  
    "EGFR":["ENSG00000146648"],      # Frequently amplified in GBM
    "IDH1":["ENSG00000138413"]
    }

    
    return h1_ensembl_map

# -------------------------------------------------
# Load MAF - improved version
# -------------------------------------------------
def load_maf():
    logging.info("Searching for GBM MAF file...")
    
    # Try different MAF file types and names
    search_terms = [
        ("Masked Somatic Mutation", None),
        ("Simple Nucleotide Variation", None),
        ("Aligned Reads", "MAF"),
        ("Simple Somatic Mutation", None)
    ]
    
    maf_files = []
    for data_type, data_format in search_terms:
        logging.info(f"Searching for data_type: {data_type}, format: {data_format}")
        files = search_gdc_files(data_type, data_format)
        
        # Filter for MAF-like files
        for f in files:
            fname = f.get('file_name', '').lower()
            if any(term in fname for term in ['.maf', 'mutation', 'somatic']):
                maf_files.append(f)
        
        if maf_files:
            break
    
    if not maf_files:
        # Try a broader search
        params = {
            "filters": {
                "op": "and",
                "content": [
                    {"op": "in", "content": {"field": "cases.project.project_id", "value": [PROJECT_ID]}},
                    {"op": "in", "content": {"field": "data_category", "value": ["Simple Nucleotide Variation"]}}
                ]
            },
            "fields": "file_id,file_name,data_format",
            "format": "JSON",
            "size": 1000
        }
        
        try:
            r = requests.post(f"{GDC_BASE}/files", json=params, timeout=30)
            r.raise_for_status()
            all_files = r.json()["data"]["hits"]
            maf_files = [f for f in all_files if '.maf' in f.get('file_name', '').lower()]
        except:
            pass
    
    if not maf_files:
        raise RuntimeError("No MAF files found for TCGA-GBM")
    
    logging.info(f"Found {len(maf_files)} potential MAF files")
    
    # Try each file until one works
    for maf_file in maf_files[:5]:
        try:
            file_name = maf_file['file_name']
            file_id = maf_file['file_id']
            
            logging.info(f"Attempting to download: {file_name}")
            content = download_file(file_id)
            
            if content is None:
                continue
            
            # Try to read the file
            try:
                if file_name.endswith('.gz'):
                    with gzip.open(BytesIO(content), 'rt', encoding='utf-8') as f:
                        # Skip comment lines and find header
                        line = f.readline()
                        while line.startswith('#'):
                            line = f.readline()
                        
                        # Now we should be at the header
                        f.seek(0)  # Go back to beginning
                        maf = pd.read_csv(f, sep='\t', comment='#', low_memory=False)
                else:
                    with BytesIO(content) as f:
                        maf = pd.read_csv(f, sep='\t', comment='#', low_memory=False)
                
                # Check for required columns
                required_cols = ['Hugo_Symbol', 'Tumor_Sample_Barcode', 'Variant_Classification']
                missing_cols = [col for col in required_cols if col not in maf.columns]
                
                if missing_cols:
                    logging.warning(f"File missing columns {missing_cols}, trying next file")
                    continue
                
                logging.info(f"Successfully loaded MAF: {file_name}")
                logging.info(f"MAF shape: {maf.shape}")
                logging.info(f"MAF columns: {list(maf.columns)}")
                logging.info(f"Sample Hugo_Symbol values: {maf['Hugo_Symbol'].dropna().unique()[:10]}")
                
                # Check if H1 genes are present
                h1_in_maf = maf[maf['Hugo_Symbol'].isin(H1_GENES)]
                if not h1_in_maf.empty:
                    logging.info(f"Found H1 genes in MAF: {h1_in_maf['Hugo_Symbol'].unique()}")
                else:
                    logging.info("No H1 genes found in this MAF file")
                
                return maf
                
            except Exception as e:
                logging.warning(f"Error reading {file_name}: {e}")
                continue
                
        except Exception as e:
            logging.warning(f"Failed to process {maf_file.get('file_name', 'unknown')}: {e}")
            continue
    
    raise RuntimeError("Could not load any valid MAF file")

# -------------------------------------------------
# Load GISTIC CNA with gene symbols
# -------------------------------------------------
def load_gistic():
    logging.info("Searching for GISTIC CNA file with gene symbols...")
    
    # Try different approaches to get CNA data with gene symbols
    
    # APPROACH 1: Look for thresholded by gene files (these usually have gene symbols)
    gistic_files = search_gdc_files("Gene Level Copy Number Scores")
    
    if not gistic_files:
        gistic_files = search_gdc_files("Copy Number Segment")
    
    # Filter for files likely to contain gene symbols
    candidate_files = []
    for f in gistic_files:
        fname = f.get('file_name', '').lower()
        # Look for files that might contain gene symbols
        if any(term in fname for term in ['all_thresholded', 'all_data', 'scores', 'by_gene']):
            candidate_files.append(f)
    
    if not candidate_files:
        candidate_files = gistic_files  # Try all files if no filtered ones
    
    logging.info(f"Found {len(candidate_files)} candidate CNA files")
    
    # Try each candidate file
    for gistic_file in candidate_files[:5]:
        try:
            file_name = gistic_file['file_name']
            file_id = gistic_file['file_id']
            
            logging.info(f"Attempting to download CNA: {file_name}")
            content = download_file(file_id)
            
            if content is None:
                continue
            
            # Try different parsing strategies
            try:
                # Try gzipped first
                if file_name.endswith('.gz'):
                    with gzip.open(BytesIO(content), 'rt') as f:
                        # Try to read as DataFrame
                        cna = pd.read_csv(f, sep='\t', low_memory=False)
                else:
                    cna = pd.read_csv(BytesIO(content), sep='\t', low_memory=False)
                
                logging.info(f"Successfully read CNA file: {file_name}")
                logging.info(f"CNA shape: {cna.shape}")
                logging.info(f"CNA columns: {list(cna.columns)}")
                
                # Try to identify gene column
                gene_col = None
                for col in cna.columns:
                    if any(term in col.lower() for term in ['gene', 'symbol', 'hugo', 'name']):
                        gene_col = col
                        break
                
                if gene_col:
                    # Set gene column as index
                    cna = cna.set_index(gene_col)
                    logging.info(f"Using '{gene_col}' as gene identifier")
                    
                    # Check if H1 genes are present
                    h1_in_cna = cna.index[cna.index.isin(H1_GENES)]
                    if len(h1_in_cna) > 0:
                        logging.info(f"Found H1 genes in CNA: {list(h1_in_cna)}")
                    else:
                        # Try to find similar names
                        logging.info(f"Gene identifiers in CNA (first 20): {cna.index[:20].tolist()}")
                        
                        # Check if index has ENSEMBL IDs and we need mapping
                        if cna.index[0].startswith('ENSG'):
                            logging.info("CNA file uses ENSEMBL IDs. Will try to map to gene symbols.")
                            # We'll handle this in the processing step
                
                return cna
                
            except Exception as e:
                logging.warning(f"Error parsing {file_name}: {e}")
                continue
                
        except Exception as e:
            logging.warning(f"Failed to process {file_name}: {e}")
            continue
    
    # If we get here, try a different approach
    logging.warning("Could not find standard CNA file, trying alternative approach...")
    
    # APPROACH 2: Try to download a specific known file
    # Sometimes TCGA-GBM has specific CNA files
    try:
        # Try to get the firehose legacy CNA data
        legacy_url = "https://gdac.broadinstitute.org/runs/analyses__2016_01_28/data/GBM/20160128/gdac.broadinstitute.org_GBM.Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3.2016012800.0.0.tar.gz"
        
        # Or use the GDC manifest approach
        manifest_params = {
            "filters": {
                "op": "and",
                "content": [
                    {"op": "in", "content": {"field": "cases.project.project_id", "value": [PROJECT_ID]}},
                    {"op": "in", "content": {"field": "data_category", "value": ["Copy Number Variation"]}},
                    {"op": "in", "content": {"field": "data_type", "value": ["Copy Number Segment", "Gene Level Copy Number"]}}
                ]
            },
            "format": "TSV",
            "fields": "file_id,file_name",
            "size": 100
        }
        
        r = requests.post(f"{GDC_BASE}/files", json=manifest_params, timeout=30)
        if r.status_code == 200:
            files_data = r.json()["data"]["hits"]
            if files_data:
                # Download the first file
                file_id = files_data[0]['file_id']
                content = download_file(file_id)
                if content:
                    cna = pd.read_csv(BytesIO(content), sep='\t')
                    return cna
    except:
        pass
    
    raise RuntimeError("Could not load CNA data")

# -------------------------------------------------
# Map ENSEMBL IDs to gene symbols for CNA data
# -------------------------------------------------
def map_ensembl_to_symbols(cna_df):
    """Try to map ENSEMBL IDs to gene symbols"""
    
    # First check if index already contains gene symbols
    if any(gene in str(idx).upper() for gene in [g.upper() for g in H1_GENES] for idx in cna_df.index[:100]):
        return cna_df
    
    # # If index looks like ENSEMBL IDs
    # if cna_df.index[0].startswith('ENSG'):
    #     logging.info("Attempting to map ENSEMBL IDs to gene symbols...")

    if str(cna_df['Chromosome'].iloc[0]).startswith('ENSG'):
        logging.info("Attempting to map ENSEMBL IDs...")
        
        # Try to get a mapping
        try:
            # Use biomart or local mapping
            mapping_url = "https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=md_ensembl_id&status=Approved&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit"
            response = requests.get(mapping_url, timeout=30)
            
            if response.status_code == 200:
                mapping_df = pd.read_csv(BytesIO(response.content), sep='\t')
                mapping_dict = {}
                for _, row in mapping_df.iterrows():
                    ensembl_id = row.get('Ensembl gene ID')
                    gene_symbol = row.get('Approved symbol')
                    if pd.notna(ensembl_id) and pd.notna(gene_symbol):
                        mapping_dict[ensembl_id] = gene_symbol
                
                # Map the CNA DataFrame
                mapped_index = []
                for idx in cna_df.index:
                    mapped_index.append(mapping_dict.get(idx, idx))
                
                cna_df.index = mapped_index
                logging.info(f"After mapping, found H1 genes: {[gene for gene in H1_GENES if gene in cna_df.index]}")
            else:
                logging.warning("Could not download gene mapping")
        except:
            logging.warning("Gene mapping failed, using ENSEMBL IDs")
    
    return cna_df

# -------------------------------------------------
# Process CNA data for H1 genes
# -------------------------------------------------
def process_cna_data(cna_df):
    """Process CNA dataframe to extract H1 gene data"""
    
    # First try to map ENSEMBL IDs to symbols
    cna_df = map_ensembl_to_symbols(cna_df)
    
    # Clean up index
    cna_df.index = cna_df.index.astype(str).str.strip()
    
    # Find H1 genes (case-insensitive)
    available_h1_genes = []
    for h1_gene in H1_GENES:
        # Exact match
        if h1_gene in cna_df.index:
            available_h1_genes.append(h1_gene)
        else:
            # Case-insensitive partial match
            matches = [idx for idx in cna_df.index if h1_gene.upper() == idx.upper()]
            if matches:
                available_h1_genes.append(matches[0])
    
    if available_h1_genes:
        logging.info(f"Found H1 genes in CNA: {available_h1_genes}")
        cna_h1 = cna_df.loc[available_h1_genes]
    else:
        logging.warning(f"No H1 genes found in CNA data. Index sample: {cna_df.index[:20].tolist()}")
        # Create empty DataFrame
        cna_h1 = pd.DataFrame(index=H1_GENES)
        if len(cna_df.columns) > 0:
            for col in cna_df.columns:
                cna_h1[col] = 0
    
    return cna_h1

# -------------------------------------------------
# Main Pipeline
# -------------------------------------------------
def run_gdc_pipeline():
    try:
        # ---- MUTATIONS ----
        maf = load_maf()
        
        # Clean up gene names
        maf['Hugo_Symbol'] = maf['Hugo_Symbol'].astype(str).str.strip()
        
        # Check what genes are actually in the MAF
        unique_genes = maf['Hugo_Symbol'].unique()
        logging.info(f"Total unique genes in MAF: {len(unique_genes)}")
        logging.info(f"Sample genes in MAF: {unique_genes[:20]}")
        
        # Filter for H1 genes
        maf_h1 = maf[maf["Hugo_Symbol"].isin(H1_GENES)]
        
        if maf_h1.empty:
            logging.warning("No H1 gene mutations found in MAF data")
            logging.info("Checking for similar gene names...")
            # Look for variations of H1 gene names
            for h1_gene in H1_GENES:
                similar = [g for g in unique_genes if h1_gene in g or h1_gene.replace('-', '') in g.replace('-', '')]
                if similar:
                    logging.info(f"Potential matches for {h1_gene}: {similar}")
            
            # Create empty mutation summary
            mutation_summary = pd.DataFrame({
                'Gene': H1_GENES,
                'Mutation_Count': 0,
                'Unique_Samples': 0,
                'Missense': 0,
                'Truncating': 0
            })
        else:
            logging.info(f"Found {len(maf_h1)} mutations in H1 genes")
            mutation_summary = (
                maf_h1.groupby("Hugo_Symbol")
                .agg(
                    Mutation_Count=("Tumor_Sample_Barcode", "count"),
                    Unique_Samples=("Tumor_Sample_Barcode", "nunique"),
                    Missense=("Variant_Classification", lambda x: (x == "Missense_Mutation").sum()),
                    Truncating=("Variant_Classification", lambda x: x.isin(
                        ["Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site"]
                    ).sum())
                )
                .reset_index()
                .rename(columns={"Hugo_Symbol": "Gene"})
            )
            
            # Add any missing H1 genes with zero values
            missing_genes = set(H1_GENES) - set(mutation_summary['Gene'])
            for gene in missing_genes:
                mutation_summary = pd.concat([
                    mutation_summary,
                    pd.DataFrame([{
                        'Gene': gene,
                        'Mutation_Count': 0,
                        'Unique_Samples': 0,
                        'Missense': 0,
                        'Truncating': 0
                    }])
                ], ignore_index=True)

        # ---- CNA ----
        cna = load_gistic()
        cna_h1 = process_cna_data(cna)
        
        logging.info(f"CNA data for H1 genes shape: {cna_h1.shape}")
        
        if cna_h1.empty or len(cna_h1.columns) == 0:
            logging.warning("No CNA data for H1 genes found")
            cna_summary = pd.DataFrame({
                'Gene': H1_GENES,
                'CNA_Altered_Samples': 0
            })
        else:
            # Create binary matrix (1 if |value| >= 1, else 0)
            # Use map() instead of applymap() to avoid FutureWarning
            cna_binary = cna_h1.copy()
            for col in cna_binary.columns:
                cna_binary[col] = cna_binary[col].map(
                    lambda x: 1 if isinstance(x, (int, float)) and abs(x) >= 1 else 0
                )
            
            cna_summary = cna_binary.sum(axis=1).reset_index()
            cna_summary.columns = ["Gene", "CNA_Altered_Samples"]
            
            # Add missing genes with zero values
            missing_genes = set(H1_GENES) - set(cna_summary['Gene'])
            for gene in missing_genes:
                cna_summary = pd.concat([
                    cna_summary,
                    pd.DataFrame([{'Gene': gene, 'CNA_Altered_Samples': 0}])
                ], ignore_index=True)

        # ---- INTEGRATION ----
        final = pd.merge(mutation_summary, cna_summary, on="Gene", how="outer").fillna(0)
        
        # Ensure all H1 genes are present
        for gene in H1_GENES:
            if gene not in final['Gene'].values:
                final = pd.concat([
                    final,
                    pd.DataFrame([{
                        'Gene': gene,
                        'Mutation_Count': 0,
                        'Unique_Samples': 0,
                        'Missense': 0,
                        'Truncating': 0,
                        'CNA_Altered_Samples': 0
                    }])
                ], ignore_index=True)
        
        # Calculate total altered samples
        final["Total_Altered_Samples"] = final["Unique_Samples"] + final["CNA_Altered_Samples"]
        final = final.sort_values("Total_Altered_Samples", ascending=False)

        # ---- SAVE ----
        final_file = os.path.join(OUTPUT_DIR, "gbm_H1_gdc_alteration_summary.csv")
        final.to_csv(final_file, index=False)

        logging.info(f"Saved GDC GBM alteration summary to {final_file}")
        logging.info("\n" + "="*50)
        logging.info("FINAL RESULTS:")
        logging.info("="*50)
        logging.info(f"\n{final.to_string()}")
        
        # Also save a more detailed version
        summary_file = os.path.join(OUTPUT_DIR, "gbm_H1_detailed_summary.txt")
        with open(summary_file, 'w') as f:
            f.write("GBM H1 Gene Alteration Summary\n")
            f.write("="*50 + "\n\n")
            f.write(final.to_string())
            f.write(f"\n\nTotal samples with any H1 alteration: {final['Total_Altered_Samples'].sum()}")
            f.write(f"\nGenes with alterations: {(final['Total_Altered_Samples'] > 0).sum()} out of {len(H1_GENES)}")
        
        logging.info(f"Detailed summary saved to: {summary_file}")
        
        return final
        
    except Exception as e:
        logging.error(f"Error in pipeline: {str(e)}", exc_info=True)
        raise

# -------------------------------------------------
# Entry Point
# -------------------------------------------------
if __name__ == "__main__":
    logging.info("Starting GDC-based GBM alteration pipeline")
    logging.info(f"Looking for H1 genes: {H1_GENES}")
    try:
        result = run_gdc_pipeline()
        logging.info("Pipeline completed successfully!")
    except Exception as e:
        logging.error(f"Pipeline failed: {e}")
        logging.error("Check if files exist in GDC and if gene names are correct.")