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

# Commonly mutated genes in Glioma/GBM (reverted to original set)
H1_GENES = [
    "TP53",      # Very common in GBM
    "PTEN",      # Common in GBM
    "EGFR",      # Frequently amplified in GBM
    "IDH1",      # Common in lower grade glioma
    "NF1",       # Common in GBM
    "PIK3CA",    # Common in many cancers
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

    # H1 gene ENSEMBL IDs for commonly mutated GBM genes
    h1_ensembl_map = {
        "TP53": ["ENSG00000141510"],
        "PTEN": ["ENSG00000171862"],
        "EGFR": ["ENSG00000146648"],
        "IDH1": ["ENSG00000138413"],
        "NF1": ["ENSG00000132155"],
        "PIK3CA": ["ENSG00000121879"]
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
    gistic_files = []

    # APPROACH 1: Prioritize 'Gene Level Copy Number Scores'
    logging.info("Attempting to find 'Gene Level Copy Number Scores'...")
    gistic_files.extend(search_gdc_files("Gene Level Copy Number Scores"))

    # APPROACH 2: Next, prioritize 'Gene Level Copy Number'
    if not gistic_files:
        logging.info("No 'Gene Level Copy Number Scores' found. Attempting to find 'Gene Level Copy Number'...")
        gistic_files.extend(search_gdc_files("Gene Level Copy Number"))

    # APPROACH 3: If still no gene-level specific files, try broader 'Copy Number Variation' category
    if not gistic_files:
        logging.info("No direct 'Gene Level Copy Number' files found. Searching broader 'Copy Number Variation' category.")
        params = {
            "filters": {
                "op": "and",
                "content": [
                    {"op": "in", "content": {"field": "cases.project.project_id", "value": [PROJECT_ID]}},
                    {"op": "in", "content": {"field": "data_category", "value": ["Copy Number Variation"]}}
                ]
            },
            "fields": "file_id,file_name,data_format",
            "format": "JSON",
            "size": 1000
        }
        try:
            r = requests.post(f"{GDC_BASE}/files", json=params, timeout=30)
            r.raise_for_status()
            all_cna_files = r.json()["data"]["hits"]
            # Filter for files likely to contain gene symbols or aggregated data
            gistic_files = [
                f for f in all_cna_files
                if any(term in f.get('file_name', '').lower() for term in ['thresholded', 'by_gene', 'scores', 'gene_expression_matrix'])
                or f.get('data_type') in ['Gene Level Copy Number', 'Gene Level Copy Number Scores'] # Ensure these are included if they appear in broader search
                or f.get('data_format') == 'TXT' # Sometimes gene-level data is just TXT
            ]
        except requests.exceptions.RequestException as e:
            logging.error(f"Error searching for Copy Number Variation files: {e}")
            gistic_files = []

    if not gistic_files:
        logging.warning("No candidate GISTIC CNA files with gene symbols found.")
        return None # Explicitly return None if no suitable files are found

    logging.info(f"Found {len(gistic_files)} candidate CNA files")

    # Try each candidate file
    for gistic_file in gistic_files[:5]: # Limit to first 5 for efficiency
        try:
            file_name = gistic_file['file_name']
            file_id = gistic_file['file_id']
            data_type = gistic_file.get('data_type')

            logging.info(f"Attempting to download CNA: {file_name} (Data Type: {data_type})")
            content = download_file(file_id)

            if content is None:
                continue

            # Try different parsing strategies
            try:
                df = None
                # Try gzipped first
                if file_name.endswith('.gz'):
                    with gzip.open(BytesIO(content), 'rt') as f:
                        df = pd.read_csv(f, sep='\t', low_memory=False)
                else:
                    df = pd.read_csv(BytesIO(content), sep='\t', low_memory=False)

                if df is None:
                    logging.warning(f"Could not read {file_name}. Skipping.")
                    continue

                logging.info(f"Successfully read CNA file: {file_name}")
                logging.info(f"CNA shape: {df.shape}")
                logging.info(f"CNA columns: {list(df.columns)}")

                # Try to identify gene column
                gene_col = None
                potential_gene_cols = ['gene', 'symbol', 'hugo', 'name', 'gene_name', 'gene_symbol', 'Gene Symbol']
                for col in df.columns:
                    if any(term in str(col).lower() for term in potential_gene_cols):
                        gene_col = col
                        break

                if gene_col:
                    # Set gene column as index
                    cna = df.set_index(gene_col)
                    # Drop any rows where the gene index is NaN or empty after setting index
                    cna = cna[cna.index.notna() & (cna.index != '')]
                    logging.info(f"Using '{gene_col}' as gene identifier. New CNA shape: {cna.shape}")

                    # Check if H1 genes are present (after setting index)
                    h1_in_cna = [gene for gene in H1_GENES if gene in cna.index]
                    if len(h1_in_cna) > 0:
                        logging.info(f"Found H1 genes in CNA data: {list(h1_in_cna)}")
                    else:
                        logging.info(f"No H1 genes found in this CNA file after gene identification. Sample index: {cna.index[:20].tolist()}")

                    return cna # Return the gene-indexed DataFrame
                else:
                    # Gracefully handle segment-level files (no direct gene column)
                    if data_type == 'Copy Number Segment' or any(term in file_name.lower() for term in ['seg', 'segment']):
                        logging.info(f"File {file_name} appears to be a segment-level file (data_type: {data_type}). It cannot be processed as gene-level directly. Skipping.")
                    else:
                        logging.warning(f"No suitable gene column found in {file_name}. Skipping this file.")
                    continue # Try next file if no gene column is identified

            except Exception as e:
                logging.warning(f"Error parsing {file_name}: {e}")
                continue

        except Exception as e:
            logging.warning(f"Failed to process {file_name}: {e}")
            continue

    logging.warning("Could not load any valid gene-level CNA file after checking multiple candidates.")
    return None # Return None if no gene-level file could be successfully loaded after all attempts.

def map_ensembl_to_symbols(cna_df):
    """Try to map ENSEMBL IDs to gene symbols"""

    # First check if index already contains gene symbols
    # Ensure index is string type for comparison
    if not cna_df.empty:
        cna_df.index = cna_df.index.astype(str)
        # Check if any H1 gene (case-insensitive) is already in the index
        if any(gene.upper() in idx.upper() for gene in H1_GENES for idx in cna_df.index):
            logging.info("Index already contains recognizable gene symbols. Skipping ENSEMBL mapping.")
            return cna_df

    # If index looks like ENSEMBL IDs (check if the first element starts with 'ENSG')
    if not cna_df.empty and isinstance(cna_df.index, pd.Index) and str(cna_df.index[0]).startswith('ENSG'):
        logging.info("Attempting to map ENSEMBL IDs in index to gene symbols...")

        try:
            mapping_url = "https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=md_ensembl_id&status=Approved&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit"
            response = requests.get(mapping_url, timeout=30, verify=False)

            if response.status_code == 200:
                raw_content = BytesIO(response.content)
                # Read the header line explicitly to inspect it
                header_line = raw_content.readline().decode('utf-8').strip()
                logging.info(f"Raw header line from mapping file: '{header_line}'")
                raw_content.seek(0) # Reset stream position

                mapping_df = pd.read_csv(raw_content, sep='\t')
                logging.info(f"Downloaded mapping file shape: {mapping_df.shape}")
                logging.info(f"Columns identified by pandas: {list(mapping_df.columns)}")
                logging.info(f"First 5 rows of mapping_df:\n{mapping_df.head().to_string()}")

                # CORRECTED COLUMN NAME based on inspection of header_line
                expected_ensembl_col = 'Ensembl gene ID'
                if expected_ensembl_col in mapping_df.columns:
                    logging.info(f"SUCCESS: Found expected Ensembl column: '{expected_ensembl_col}'")
                    non_null_ensembl_ids = mapping_df[expected_ensembl_col].dropna()
                    if not non_null_ensembl_ids.empty:
                        logging.info(f"Sample of non-null Ensembl IDs from mapping file: {non_null_ensembl_ids.head().tolist()}")
                    else:
                        logging.warning(f"Column '{expected_ensembl_col}' found, but contains no non-null Ensembl IDs.")
                else:
                    logging.error(f"CRITICAL ERROR: Column '{expected_ensembl_col}' NOT found in downloaded mapping file. Available columns: {mapping_df.columns.tolist()}")

                mapping_dict = {}
                for _, row in mapping_df.iterrows():
                    ensembl_id = row.get(expected_ensembl_col)
                    gene_symbol = row.get('Approved symbol')
                    if pd.notna(ensembl_id) and pd.notna(gene_symbol):
                        # Ensure Ensembl ID is stripped of version number for robust mapping
                        mapping_dict[str(ensembl_id).split('.')[0].strip()] = gene_symbol.strip()
                logging.info(f"Size of mapping dictionary: {len(mapping_dict)}")
                logging.info(f"Sample keys from mapping_dict: {list(mapping_dict.keys())[:10]}")

                h1_ensembl_map_local = get_gene_symbol_mapping()
                for h1_gene_symbol in H1_GENES:
                    ensembl_ids_for_h1 = h1_ensembl_map_local.get(h1_gene_symbol, [])
                    for ensembl_id in ensembl_ids_for_h1:
                        # Strip version for lookup in our map as well
                        clean_ensembl_id = ensembl_id.split('.')[0].strip() if ensembl_id else None
                        if clean_ensembl_id and clean_ensembl_id in mapping_dict:
                            logging.info(f"Found '{clean_ensembl_id}' for {h1_gene_symbol} in mapping_dict. Maps to: {mapping_dict[clean_ensembl_id]}")
                        else:
                            logging.warning(f"Our ENSEMBL {ensembl_id} for {h1_gene_symbol} NOT found in mapping_dict (or is None). Check if correct mapping or if it has a version number.")

                mapped_index = []
                for idx in cna_df.index:
                    # Strip any version numbers from ENSEMBL IDs like 'ENSG00000223972.5' -> 'ENSG00000223972'
                    clean_idx = idx.split('.')[0].strip()
                    mapped_symbol = mapping_dict.get(clean_idx, idx) # Map using cleaned ID, fallback to original if not found
                    mapped_index.append(mapped_symbol)

                cna_df.index = mapped_index
                logging.info(f"After mapping, found H1 genes: {[gene for gene in H1_GENES if gene in cna_df.index]}")
            else:
                logging.warning(f"Could not download gene mapping. Status code: {response.status_code}")
        except requests.exceptions.RequestException as e:
            logging.warning(f"Gene mapping failed due to request error: {e}")
        except Exception as e:
            logging.warning(f"Gene mapping failed due to an unexpected error: {e}", exc_info=True)

    return cna_df


# -------------------------------------------------
# Process CNA data for H1 genes
# -------------------------------------------------
def process_cna_data(cna_df):
    """Process CNA dataframe to extract H1 gene data"""

    if cna_df is None or cna_df.empty:
        logging.warning("No CNA data provided or data is empty. Returning empty DataFrame for H1 genes.")
        cna_h1 = pd.DataFrame(index=H1_GENES)
        # Add a dummy numerical column if needed for later processing (e.g., sum)
        if H1_GENES: # Ensure H1_GENES is not empty
            cna_h1['CNA_Score'] = 0
        return cna_h1

    # Make a copy to avoid SettingWithCopyWarning
    processed_df = cna_df.copy()

    # Prioritize 'gene_name' column if available and not already the index
    if 'gene_name' in processed_df.columns and not processed_df.index.name == 'gene_name':
        logging.info("Setting 'gene_name' column as index.")
        processed_df = processed_df.set_index('gene_name')
        # Drop any rows where the gene index is NaN or empty after setting index
        processed_df = processed_df[processed_df.index.notna() & (processed_df.index != '')]
    
    # Now, try to map ENSEMBL IDs to symbols if the index appears to be ENSEMBL IDs
    processed_df = map_ensembl_to_symbols(processed_df)

    # Clean up index (make sure it's string and strip whitespace)
    processed_df.index = processed_df.index.astype(str).str.strip()

    # Filter for H1 genes (case-insensitive for robustness)
    # Create a mapping from cleaned H1_GENES to actual gene names in index
    h1_gene_map = {}
    for h1_gene in H1_GENES:
        # Check for exact match first
        if h1_gene in processed_df.index:
            h1_gene_map[h1_gene] = h1_gene
        else:
            # Check for case-insensitive match
            matches = [idx for idx in processed_df.index if idx.lower() == h1_gene.lower()]
            if matches:
                h1_gene_map[h1_gene] = matches[0]

    found_h1_gene_symbols = list(h1_gene_map.values())

    if found_h1_gene_symbols:
        logging.info(f"Found H1 genes in CNA: {list(h1_gene_map.keys())}")
        cna_h1 = processed_df.loc[found_h1_gene_symbols]
        # Reindex to ensure all H1_GENES are present in the final output, filling missing with 0
        # This also ensures the order matches H1_GENES
        cna_h1 = cna_h1.reindex(list(h1_gene_map.keys()), fill_value=0.0)
    else:
        logging.warning(f"No H1 genes found in CNA data after processing. Index sample: {processed_df.index[:20].tolist()}")
        # Create empty DataFrame with H1_GENES as index
        cna_h1 = pd.DataFrame(index=H1_GENES)
        
    # Drop any non-numeric columns if they are not scores (e.g., Chromosome, Start, End)
    numeric_cols = cna_h1.select_dtypes(include=np.number).columns
    if not numeric_cols.empty:
        cna_h1 = cna_h1[numeric_cols]
    else:
        # If no numeric columns found, create one with zeros
        logging.warning("No numeric columns found in cna_h1, creating 'CNA_Score' column with zeros.")
        cna_h1['CNA_Score'] = 0.0

    return cna_h1


# -------------------------------------------------
# Main Pipeline
# -------------------------------------------------
def run_gdc_pipeline():
    try:
        # ---- MUTATIONS ----
        maf = load_maf()

        maf['Hugo_Symbol'] = maf['Hugo_Symbol'].astype(str).str.strip()

        unique_genes = maf['Hugo_Symbol'].unique()
        logging.info(f"Total unique genes in MAF: {len(unique_genes)}")
        logging.info(f"Sample genes in MAF: {unique_genes[:20]}")

        maf_h1 = maf[maf["Hugo_Symbol"].isin(H1_GENES)]

        if maf_h1.empty:
            logging.warning("No H1 gene mutations found in MAF data")
            logging.info("Checking for similar gene names...")
            for h1_gene in H1_GENES:
                similar = [g for g in unique_genes if h1_gene in g or h1_gene.replace('-', '') in g.replace('-', '')]
                if similar:
                    logging.info(f"Potential matches for {h1_gene}: {similar}")

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
            cna_binary = cna_h1.copy()
            for col in cna_binary.columns:
                cna_binary[col] = cna_binary[col].map(
                    lambda x: 1 if isinstance(x, (int, float)) and abs(x) >= 1 else 0
                )

            cna_summary = cna_binary.sum(axis=1).reset_index()
            cna_summary.columns = ["Gene", "CNA_Altered_Samples"]

            missing_genes = set(H1_GENES) - set(cna_summary['Gene'])
            for gene in missing_genes:
                cna_summary = pd.concat([
                    cna_summary,
                    pd.DataFrame([{'Gene': gene, 'CNA_Altered_Samples': 0}])
                ], ignore_index=True)

        # ---- INTEGRATION ----
        final = pd.merge(mutation_summary, cna_summary, on="Gene", how="outer").fillna(0)

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

        final["Total_Altered_Samples"] = final["Unique_Samples"] + final["CNA_Altered_Samples"]
        final = final.sort_values("Total_Altered_Samples", ascending=False)

        final_file = os.path.join(OUTPUT_DIR, "gbm_H1_gdc_alteration_summary.csv")
        final.to_csv(final_file, index=False)

        logging.info(f"Saved GDC GBM alteration summary to {final_file}")
        logging.info("\n" + "="*50)
        logging.info("FINAL RESULTS:")
        logging.info("="*50)
        logging.info(f"\n{final.to_string()}")

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