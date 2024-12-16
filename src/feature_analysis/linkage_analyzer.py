"""
Linkage Disequilibrium Analysis Module for identifying trait associations.
"""
import logging
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from scipy import stats

logger = logging.getLogger(__name__)

class LinkageAnalyzer:
    """Analyzes linkage disequilibrium between genetic markers."""

    def __init__(self):
        """Initialize LinkageAnalyzer."""
        self.ld_results = pd.DataFrame()

    def calculate_ld(self, genotypes: pd.DataFrame) -> pd.DataFrame:
        """
        Calculate linkage disequilibrium between markers.

        Args:
            genotypes: DataFrame with marker genotypes (rows=samples, cols=markers)

        Returns:
            DataFrame with pairwise LD statistics
        """
        try:
            n_markers = genotypes.shape[1]
            ld_data = []

            for i in range(n_markers):
                for j in range(i + 1, n_markers):
                    # Calculate allele frequencies
                    p1 = np.clip(genotypes.iloc[:, i].mean() / 2, 0, 1)
                    p2 = np.clip(genotypes.iloc[:, j].mean() / 2, 0, 1)
                    q1 = 1 - p1
                    q2 = 1 - p2

                    # Calculate haplotype frequency
                    hap_freq = np.clip(np.mean(genotypes.iloc[:, i] * genotypes.iloc[:, j]) / 4, 0, 1)

                    # Calculate D (LD coefficient)
                    d = hap_freq - (p1 * p2)

                    # Calculate D' (normalized LD)
                    d_max = min(p1 * q2, q1 * p2) if d > 0 else min(p1 * p2, q1 * q2)
                    d_prime = abs(d / d_max) if d_max != 0 else 0

                    # Calculate r² (correlation coefficient)
                    denom = p1 * q1 * p2 * q2
                    r2 = (d * d) / denom if denom != 0 else 0

                    ld_data.append({
                        'marker1': genotypes.columns[i],
                        'marker2': genotypes.columns[j],
                        'D': abs(d),  # Use absolute D value
                        'D_prime': d_prime,
                        'r2': r2
                    })

            self.ld_results = pd.DataFrame(ld_data)
            logger.info(f"Calculated LD statistics for {n_markers} markers")
            return self.ld_results

        except Exception as e:
            logger.error(f"Error calculating LD: {str(e)}")
            raise

    def identify_trait_associations(self, ld_threshold: float = 0.8) -> List[Dict]:
        """
        Identify significant trait associations based on LD.

        Args:
            ld_threshold: Threshold for significant LD (D' value)

        Returns:
            List of significant marker associations
        """
        try:
            significant = self.ld_results[self.ld_results['D_prime'] >= ld_threshold]
            associations = significant.to_dict('records')

            logger.info(f"Found {len(associations)} significant trait associations")
            return associations

        except Exception as e:
            logger.error(f"Error identifying trait associations: {str(e)}")
            raise

    def generate_ld_report(self) -> Dict:
        """
        Generate summary statistics for LD analysis.

        Returns:
            Dictionary containing LD summary statistics
        """
        try:
            stats = {
                'total_pairs': len(self.ld_results),
                'mean_d_prime': float(np.clip(self.ld_results['D_prime'].mean(), 0, 1)),
                'mean_r2': float(np.clip(self.ld_results['r2'].mean(), 0, 1)),
                'high_ld_pairs': len(self.ld_results[self.ld_results['D_prime'] >= 0.8])
            }

            logger.info("Generated LD analysis report")
            return stats

        except Exception as e:
            logger.error(f"Error generating LD report: {str(e)}")
            raise
