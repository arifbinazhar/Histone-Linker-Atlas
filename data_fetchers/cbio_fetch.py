"""
cBioPortal Mutation Fetcher for Linker Histone Analysis in Glioblastoma

Author: HistoneLinker-Atlas Team
Date: 2024
Purpose: Fetch somatic mutations in H1 family genes from TCGA GBM using cBioPortal API
         to identify recurrent mutations for chromatin biology studies.

Biological Context:
- Linker histones (H1 family) regulate higher-order chromatin structure
- Somatic mutations in H1 genes may disrupt chromatin organization in GBM
- Recurrent mutations may indicate driver events in cancer progression

Engineering Standards:
- Fault-tolerant with retry logic
- Batching for API limits (100 samples per request)
- Comprehensive logging
- Data validation and sanity checks
"""

import os
import sys
import json
import time
import logging
from typing import List, Dict, Any, Optional, Tuple
from pathlib import Path
import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import yaml

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('data_fetchers/cbio_fetch.py'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

class CBioPortalClient:
    """Robust cBioPortal API client with error handling and retry logic"""
    
    BASE_URL = "https://www.cbioportal.org/api"
    
    def __init__(self, timeout: int = 30, max_retries: int = 3):
        """Initialize API client with retry strategy"""
        self.session = self._create_session(max_retries)
        self.timeout = timeout
        self.max_retries = max_retries
        
    def _create_session(self, max_retries: int) -> requests.Session:
        """Create HTTP session with retry logic"""
        session = requests.Session()
        retry_strategy = Retry(
            total=max_retries,
            backoff_factor=1,
            status_forcelist=[429, 500, 502, 503, 504],
            allowed_methods=["GET", "POST"]
        )
        adapter = HTTPAdapter(max_retries=retry_strategy)
        session.mount("https://", adapter)
        return session
    
    def _make_request(self, method: str, endpoint: str, **kwargs) -> Optional[Dict]:
        """Make API request with error handling"""
        url = f"{self.BASE_URL}{endpoint}"
        
        try:
            if method.upper() == "GET":
                response = self.session.get(url, timeout=self.timeout, **kwargs)
            elif method.upper() == "POST":
                response = self.session.post(url, timeout=self.timeout, **kwargs)
            else:
                raise ValueError(f"Unsupported method: {method}")
            
            response.raise_for_status()
            
            if response.status_code == 204:
                return None
            
            return response.json()
            
        except requests.exceptions.RequestException as e:
            logger.error(f"API request failed: {e} - URL: {url}")
            raise
        except json.JSONDecodeError as e:
            logger.error(f"JSON decode error: {e} - Response: {response.text[:200]}")
            raise

class H1MutationAnalyzer:
    """Main analyzer for H1 family mutations in GBM"""
    
    # H1 family genes with aliases resolved
    H1_GENES = [
        "H1F0",      # H1.0
        "HIST1H1A",  # H1.1
        "HIST1H1B",  # H1.5
        "HIST1H1C",  # H1.2
        "HIST1H1D",  # H1.3
        "HIST1H1E",  # H1.4
        "H1FX"       # H1.10
    ]
    
    STUDY_ID = "gbm_tcga"
    
    def __init__(self, output_dir: Path, batch_size: int = 100):
        """
        Initialize mutation analyzer
        
        Args:
            output_dir: Path to processed data directory
            batch_size: Number of samples per API request (cBioPortal limit: 100-200)
        """
        self.output_dir = output_dir
        self.batch_size = min(batch_size, 200)  # Enforce API limit
        self.client = CBioPortalClient()
        self.molecular_profile_id = None
        
        # Create output directory
        self.gbm_dir = output_dir / "gbm"
        self.gbm_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"Initialized H1MutationAnalyzer for {self.STUDY_ID}")
        logger.info(f"Output directory: {self.gbm_dir}")
    
    def resolve_molecular_profile(self) -> str:
        """
        Identify the correct molecular profile for mutations in TCGA GBM
        
        Returns:
            molecularProfileId for mutation data
        
        Note: Must select profile with 'mutations' in ID and matching study
        """
        logger.info("Resolving molecular profile for mutation data...")
        
        try:
            profiles = self.client._make_request(
                "GET", 
                f"/api/molecular-profiles?studyId={self.STUDY_ID}"
            )
            
            if not profiles:
                raise ValueError(f"No molecular profiles found for study {self.STUDY_ID}")
            
            # Filter for mutation profiles specific to GBM study
            mutation_profiles = [
                p for p in profiles 
                if "mutations" in p["molecularProfileId"].lower() 
                and self.STUDY_ID in p["molecularProfileId"]
            ]
            
            if not mutation_profiles:
                raise ValueError(f"No mutation profiles found for study {self.STUDY_ID}")
            
            # Select the most appropriate profile (prioritize MC3 if available)
            selected_profile = None
            for profile in mutation_profiles:
                profile_id = profile["molecularProfileId"]
                if "mutations" in profile_id and self.STUDY_ID in profile_id:
                    # Prefer harmonized MC3 calls
                    if "mc3" in profile_id.lower():
                        selected_profile = profile
                        break
                    elif not selected_profile:  # Fallback to first mutation profile
                        selected_profile = profile
            
            self.molecular_profile_id = selected_profile["molecularProfileId"]
            logger.info(f"Selected molecular profile: {self.molecular_profile_id}")
            logger.info(f"Profile name: {selected_profile.get('name', 'N/A')}")
            
            return self.molecular_profile_id
            
        except Exception as e:
            logger.error(f"Failed to resolve molecular profile: {e}")
            raise
    
    def fetch_sample_ids(self) -> List[str]:
        """
        Fetch all sample IDs for TCGA GBM study
        
        Returns:
            List of sample IDs
        """
        logger.info(f"Fetching sample IDs for study {self.STUDY_ID}...")
        
        try:
            samples = self.client._make_request(
                "GET", 
                f"/api/studies/{self.STUDY_ID}/samples"
            )
            
            sample_ids = [sample["sampleId"] for sample in samples]
            logger.info(f"Found {len(sample_ids)} samples")
            
            # Validate sample format (should be TCGA IDs)
            tcga_samples = [sid for sid in sample_ids if sid.startswith("TCGA-")]
            if len(tcga_samples) != len(sample_ids):
                logger.warning(f"Found {len(sample_ids) - len(tcga_samples)} non-TCGA samples")
            
            return sample_ids
            
        except Exception as e:
            logger.error(f"Failed to fetch sample IDs: {e}")
            raise
    
    def fetch_mutations_batch(self, sample_batch: List[str]) -> List[Dict]:
        """
        Fetch mutations for a batch of samples
        
        Args:
            sample_batch: List of sample IDs
            
        Returns:
            List of mutation records
        """
        if not self.molecular_profile_id:
            raise ValueError("Molecular profile ID not resolved")
        
        payload = {
            "geneSymbols": self.H1_GENES,
            "sampleIds": sample_batch
        }
        
        endpoint = f"/api/molecular-profiles/{self.molecular_profile_id}/mutations/fetch"
        
        try:
            logger.info(f"Fetching mutations for {len(sample_batch)} samples...")
            mutations = self.client._make_request("POST", endpoint, json=payload)
            
            if mutations:
                logger.info(f"Retrieved {len(mutations)} mutation records")
            else:
                logger.info("No mutations found in this batch")
            
            return mutations or []
            
        except Exception as e:
            logger.error(f"Failed to fetch mutations for batch: {e}")
            # Return empty list to allow continuation
            return []
    
    def batch_samples(self, sample_ids: List[str]) -> List[List[str]]:
        """
        Split samples into batches for API requests
        
        Args:
            sample_ids: List of all sample IDs
            
        Returns:
            List of sample batches
        """
        batches = []
        for i in range(0, len(sample_ids), self.batch_size):
            batch = sample_ids[i:i + self.batch_size]
            batches.append(batch)
        
        logger.info(f"Split {len(sample_ids)} samples into {len(batches)} batches "
                   f"(batch size: {self.batch_size})")
        return batches
    
    def fetch_all_mutations(self, sample_ids: List[str]) -> pd.DataFrame:
        """
        Fetch mutations for all samples with batching
        
        Args:
            sample_ids: List of all sample IDs
            
        Returns:
            DataFrame with all mutation records
        """
        all_mutations = []
        batches = self.batch_samples(sample_ids)
        
        for i, batch in enumerate(batches, 1):
            logger.info(f"Processing batch {i}/{len(batches)} ({len(batch)} samples)")
            
            try:
                mutations = self.fetch_mutations_batch(batch)
                all_mutations.extend(mutations)
                
                # Be respectful to API - add delay between batches
                if i < len(batches):
                    time.sleep(1)
                    
            except Exception as e:
                logger.error(f"Batch {i} failed: {e}")
                # Continue with next batch
                continue
        
        if not all_mutations:
            logger.warning("No mutations found for any H1 gene")
            return pd.DataFrame()
        
        # Convert to DataFrame
        df = pd.DataFrame(all_mutations)
        logger.info(f"Total mutations collected: {len(df)}")
        
        return df
    
    def classify_mutation_type(self, mutation: Dict) -> str:
        """
        Classify mutation by functional impact
        
        Args:
            mutation: Mutation record from cBioPortal
            
        Returns:
            Mutation type classification
        """
        mutation_type = mutation.get("mutationType", "").upper()
        protein_change = mutation.get("proteinChange", "").upper()
        
        # Truncating mutations
        if any(term in mutation_type for term in ["NONSENSE", "FRAME_SHIFT", "SPLICE_SITE"]):
            return "truncating"
        
        # Missense mutations
        elif "MISSENSE" in mutation_type:
            return "missense"
        
        # In-frame mutations
        elif "IN_FRAME" in mutation_type:
            return "in_frame"
        
        # Silent/synonymous
        elif any(term in mutation_type for term in ["SILENT", "SYNONYMOUS"]):
            return "silent"
        
        # Non-coding
        elif any(term in mutation_type for term in ["NONSTOP", "NONSTART", "5UTR", "3UTR"]):
            return "non_coding"
        
        else:
            return "other"
    
    def summarize_mutations(self, mutations_df: pd.DataFrame) -> pd.DataFrame:
        """
        Generate summary statistics for H1 mutations
        
        Args:
            mutations_df: Raw mutation DataFrame
            
        Returns:
            Summary DataFrame with gene-level statistics
        """
        if mutations_df.empty:
            logger.warning("No mutations to summarize")
            return pd.DataFrame()
        
        summary_data = []
        
        for gene in self.H1_GENES:
            gene_mutations = mutations_df[mutations_df["geneSymbol"] == gene]
            
            # Basic counts
            mutation_count = len(gene_mutations)
            unique_samples = gene_mutations["sampleId"].nunique() if mutation_count > 0 else 0
            
            # Mutation type classification
            mutation_types = []
            for _, mut in gene_mutations.iterrows():
                mut_type = self.classify_mutation_type(mut.to_dict())
                mutation_types.append(mut_type)
            
            # Count by type
            type_counts = pd.Series(mutation_types).value_counts()
            
            summary_data.append({
                "Gene": gene,
                "Mutation_Count": mutation_count,
                "Unique_Samples": unique_samples,
                "Missense_Count": type_counts.get("missense", 0),
                "Truncating_Count": type_counts.get("truncating", 0),
                "InFrame_Count": type_counts.get("in_frame", 0),
                "Silent_Count": type_counts.get("silent", 0),
                "NonCoding_Count": type_counts.get("non_coding", 0),
                "Other_Count": type_counts.get("other", 0)
            })
        
        summary_df = pd.DataFrame(summary_data)
        
        # Calculate mutation frequency
        total_samples = mutations_df["sampleId"].nunique()
        if total_samples > 0:
            summary_df["Mutation_Frequency"] = summary_df["Unique_Samples"] / total_samples
            summary_df["Mutation_Frequency_Percent"] = summary_df["Mutation_Frequency"] * 100
        
        # Sort by mutation count
        summary_df = summary_df.sort_values("Mutation_Count", ascending=False)
        
        return summary_df
    
    def save_results(self, raw_mutations: pd.DataFrame, summary: pd.DataFrame):
        """
        Save mutation data to files
        
        Args:
            raw_mutations: Raw mutation DataFrame
            summary: Summary DataFrame
        """
        # Save raw mutations
        raw_path = self.gbm_dir / "gbm_H1_mutations_raw.csv"
        if not raw_mutations.empty:
            raw_mutations.to_csv(raw_path, index=False)
            logger.info(f"Saved raw mutations to {raw_path}")
            
            # Log sample mutation statistics
            mutated_samples = raw_mutations["sampleId"].nunique()
            total_samples = len(set().union(*[batch for batch in self.batch_samples(
                list(raw_mutations["sampleId"].unique())
            )]))
            logger.info(f"Mutations found in {mutated_samples}/{total_samples} samples")
        
        # Save summary
        summary_path = self.gbm_dir / "gbm_H1_mutation_summary.csv"
        if not summary.empty:
            summary.to_csv(summary_path, index=False)
            logger.info(f"Saved mutation summary to {summary_path}")
            
            # Log top mutated genes
            top_genes = summary.head(3)
            for _, row in top_genes.iterrows():
                logger.info(f"{row['Gene']}: {row['Mutation_Count']} mutations "
                          f"({row.get('Unique_Samples', 'N/A')} samples)")
    
    def run(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Execute complete mutation analysis pipeline
        
        Returns:
            Tuple of (raw_mutations_df, summary_df)
        """
        logger.info("=" * 60)
        logger.info("Starting H1 Mutation Analysis Pipeline")
        logger.info("=" * 60)
        
        try:
            # Step 1: Resolve molecular profile
            profile_id = self.resolve_molecular_profile()
            
            # Step 2: Fetch all sample IDs
            sample_ids = self.fetch_sample_ids()
            
            if not sample_ids:
                raise ValueError("No sample IDs found")
            
            # Step 3: Fetch mutations with batching
            raw_mutations = self.fetch_all_mutations(sample_ids)
            
            # Step 4: Generate summary
            summary = self.summarize_mutations(raw_mutations)
            
            # Step 5: Save results
            self.save_results(raw_mutations, summary)
            
            logger.info("=" * 60)
            logger.info("Mutation analysis completed successfully")
            logger.info("=" * 60)
            
            return raw_mutations, summary
            
        except Exception as e:
            logger.error(f"Pipeline failed: {e}")
            raise

def main():
    """Main entry point for the module"""
    # Read configuration
    config_path = Path("config.yaml")
    if config_path.exists():
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
        output_dir = Path(config.get("output_dir", "data/processed"))
        batch_size = config.get("batch_size", 100)
    else:
        output_dir = Path("data/processed")
        batch_size = 100
    
    # Ensure log directory exists
    log_dir = Path("logs")
    log_dir.mkdir(exist_ok=True)
    
    try:
        analyzer = H1MutationAnalyzer(output_dir, batch_size)
        raw_mutations, summary = analyzer.run()
        
        # Print summary to console for quick review
        if not summary.empty:
            print("\n" + "=" * 60)
            print("H1 MUTATION SUMMARY (TCGA GBM)")
            print("=" * 60)
            print(summary[["Gene", "Mutation_Count", "Unique_Samples", 
                          "Missense_Count", "Truncating_Count"]].to_string(index=False))
            print("\nNote: Use for experimental prioritization of chromatin regulators")
            
    except Exception as e:
        logger.error(f"Fatal error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()2026-02-01 02:09:06,255 - __main__ - INFO - Initialized H1MutationAnalyzer for gbm_tcga
2026-02-01 02:09:06,255 - __main__ - INFO - Output directory: data\processed\gbm
2026-02-01 02:09:06,256 - __main__ - INFO - ============================================================
2026-02-01 02:09:06,256 - __main__ - INFO - Starting H1 Mutation Analysis Pipeline
2026-02-01 02:09:06,256 - __main__ - INFO - ============================================================
2026-02-01 02:09:06,257 - __main__ - INFO - Resolving molecular profile for mutation data...
2026-02-01 02:09:07,495 - __main__ - ERROR - API request failed: 404 Client Error:  for url: https://www.cbioportal.org/api/api/molecular-profiles?studyId=gbm_tcga - URL: https://www.cbioportal.org/api/api/molecular-profiles?studyId=gbm_tcga
2026-02-01 02:09:07,496 - __main__ - ERROR - Failed to resolve molecular profile: 404 Client Error:  for url: https://www.cbioportal.org/api/api/molecular-profiles?studyId=gbm_tcga
2026-02-01 02:09:07,496 - __main__ - ERROR - Pipeline failed: 404 Client Error:  for url: https://www.cbioportal.org/api/api/molecular-profiles?studyId=gbm_tcga
2026-02-01 02:09:07,497 - __main__ - ERROR - Fatal error: 404 Client Error:  for url: https://www.cbioportal.org/api/api/molecular-profiles?studyId=gbm_tcga
