"""
VCF file handler for processing variant data.
"""
from pathlib import Path
from typing import Dict, List, Optional
import vcf
import logging
import pandas as pd

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class VcfHandler:
    """Handles loading and processing of VCF files containing variant data."""

    def __init__(self, min_quality: float = 30.0):
        """Initialize the VCF handler.

        Args:
            min_quality: Minimum quality score for variants
        """
        self.min_quality = min_quality

    def load_vcf(self, filepath: Path) -> vcf.Reader:
        """Load and validate VCF file.

        Args:
            filepath: Path to the VCF file

        Returns:
            VCF Reader object

        Raises:
            FileNotFoundError: If the file doesn't exist
            ValueError: If the file is invalid
        """
        try:
            if not filepath.exists():
                raise FileNotFoundError(f"VCF file not found: {filepath}")


            logger.info(f"Loading VCF file: {filepath}")
            return vcf.Reader(open(filepath, 'r'))

        except Exception as e:
            logger.error(f"Error loading VCF file {filepath}: {str(e)}")
            raise

    def filter_variants(self, vcf_reader: vcf.Reader,
                       region: Optional[str] = None) -> pd.DataFrame:
        """Filter variants based on quality and region.

        Args:
            vcf_reader: VCF Reader object
            region: Optional region filter (e.g., "chr1:1000-2000")

        Returns:
            DataFrame containing filtered variants
        """
        variants = []

        try:
            for record in vcf_reader:
                # Quality filter
                if record.QUAL and record.QUAL < self.min_quality:
                    logger.debug(f"Filtering low quality variant: {record.ID}")
                    continue

                # Region filter
                if region:
                    chrom, pos = region.split(':')[0], int(region.split(':')[1])
                    if record.CHROM != chrom or record.POS != pos:
                        continue

                variant_data = {
                    'CHROM': record.CHROM,
                    'POS': record.POS,
                    'ID': record.ID,
                    'REF': record.REF,
                    'ALT': ','.join(str(a) for a in record.ALT),
                    'QUAL': record.QUAL,
                    'FILTER': ','.join(record.FILTER) if record.FILTER else 'PASS'
                }

                # Add INFO fields
                for field in record.INFO:
                    variant_data[f'INFO_{field}'] = record.INFO[field]

                variants.append(variant_data)

        except Exception as e:
            logger.error(f"Error processing variants: {str(e)}")
            raise

        return pd.DataFrame(variants)

    def organize_by_gene(self, variants: pd.DataFrame) -> Dict[str, pd.DataFrame]:
        """Organize variants by gene.

        Args:
            variants: DataFrame of variants

        Returns:
            Dictionary mapping gene names to variant DataFrames
        """
        gene_variants = {}

        if 'INFO_Gene' not in variants.columns:
            logger.warning("No gene information found in variants")
            return gene_variants

        for gene in variants['INFO_Gene'].unique():
            if pd.isna(gene):
                continue
            gene_variants[gene] = variants[variants['INFO_Gene'] == gene]

        return gene_variants
