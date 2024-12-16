"""
VCF file handling functionality.
"""
from pathlib import Path
from typing import List, Dict
import pandas as pd
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class VCFHandler:
    """Handles VCF file operations for variant data."""

    def __init__(self):
        """Initialize VCFHandler."""
        self.logger = logger

    def load_vcf(self, file_path: str) -> List[Dict]:
        """
        Load variants from a VCF file.

        Args:
            file_path: Path to the VCF file

        Returns:
            List of variant dictionaries
        """
        try:
            self.logger.info(f"Loading VCF file: {file_path}")
            variants = []
            with open(file_path, 'r') as vcf_file:
                for line in vcf_file:
                    if line.startswith('#'):
                        continue
                    fields = line.strip().split('\t')
                    variant = {
                        'CHROM': fields[0],
                        'POS': int(fields[1]),
                        'ID': fields[2],
                        'REF': fields[3],
                        'ALT': fields[4],
                        'QUAL': float(fields[5]) if fields[5] != '.' else 0,
                        'FILTER': fields[6],
                        'INFO': fields[7],
                        'FORMAT': fields[8],
                        'samples': fields[9:]
                    }
                    variants.append(variant)
            self.logger.info(f"Successfully loaded {len(variants)} variants")
            return variants
        except FileNotFoundError:
            self.logger.error(f"VCF file not found: {file_path}")
            raise
        except Exception as e:
            self.logger.error(f"Error loading VCF file: {str(e)}")
            raise

    def filter_variants(self, variants: List[Dict], min_qual: float = 0) -> List[Dict]:
        """
        Filter variants based on quality score.

        Args:
            variants: List of variant dictionaries
            min_qual: Minimum quality score to keep

        Returns:
            Filtered list of variant dictionaries
        """
        try:
            filtered = [var for var in variants if var['QUAL'] >= min_qual]
            self.logger.info(f"Filtered {len(variants) - len(filtered)} variants below quality {min_qual}")
            return filtered
        except Exception as e:
            self.logger.error(f"Error filtering variants: {str(e)}")
            raise

    def extract_genotypes(self, variants: List[Dict]) -> pd.DataFrame:
        """
        Extract genotype information from variants.

        Args:
            variants: List of variant dictionaries

        Returns:
            DataFrame with genotype information
        """
        try:
            genotypes = []
            for var in variants:
                # Convert genotypes to numeric format (0=ref/ref, 1=ref/alt, 2=alt/alt)
                gt_values = []
                for sample in var['samples']:
                    gt = sample.split(':')[0]  # Get GT field
                    if gt == '0/0':
                        gt_values.append(0)
                    elif gt in ['0/1', '1/0']:
                        gt_values.append(1)
                    elif gt == '1/1':
                        gt_values.append(2)
                    else:
                        gt_values.append(None)
                genotypes.append(gt_values)

            # Create DataFrame with variant IDs as index
            df = pd.DataFrame(genotypes,
                            columns=[f'SAMPLE{i+1}' for i in range(len(variants[0]['samples']))],
                            index=[var['ID'] for var in variants])

            self.logger.info(f"Extracted genotypes for {len(variants)} variants")
            return df
        except Exception as e:
            self.logger.error(f"Error extracting genotypes: {str(e)}")
            raise
