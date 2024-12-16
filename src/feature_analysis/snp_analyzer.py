"""
SNP Analysis Module for evaluating trait-associated variants.
"""
import logging
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from scipy import stats

logger = logging.getLogger(__name__)

class SNPAnalyzer:
    """Analyzes SNPs for trait associations and population differentiation."""

    def __init__(self):
        """Initialize SNPAnalyzer."""
        self.snp_stats = pd.DataFrame()

    def evaluate_trait_snps(self, genotypes: pd.DataFrame,
                          phenotypes: pd.DataFrame) -> pd.DataFrame:
        """
        Evaluate SNPs for trait associations.

        Args:
            genotypes: DataFrame with SNP genotypes
            phenotypes: DataFrame with trait phenotypes

        Returns:
            DataFrame with SNP association statistics
        """
        try:
            snp_data = []

            for snp in genotypes.columns:
                # Calculate allele frequencies by trait value
                freq_by_trait = {}
                for trait_value in phenotypes['trait_value'].unique():
                    trait_mask = phenotypes['trait_value'] == trait_value
                    freq_by_trait[trait_value] = genotypes.loc[trait_mask, snp].mean() / 2

                # Create 2x2 contingency table for each genotype
                # Combine heterozygous and homozygous alternate for binary comparison
                binary_genotypes = (genotypes[snp] > 0).astype(int)
                contingency = pd.crosstab(binary_genotypes,
                                        phenotypes['trait_value'])

                # Ensure 2x2 shape for Fisher's exact test
                if contingency.shape == (2, 2):
                    odds_ratio, pval = stats.fisher_exact(contingency)
                else:
                    odds_ratio, pval = np.nan, 1.0

                # Chi-square test
                chi2, chi2_pval = stats.chi2_contingency(contingency)[:2]

                snp_data.append({
                    'snp': snp,
                    'chi2_statistic': chi2,
                    'chi2_p_value': chi2_pval,
                    'fisher_p_value': pval,
                    'odds_ratio': odds_ratio,
                    'allele_freqs': freq_by_trait
                })

            self.snp_stats = pd.DataFrame(snp_data)
            logger.info(f"Evaluated {len(snp_data)} SNPs for trait associations")
            return self.snp_stats

        except Exception as e:
            logger.error(f"Error evaluating trait SNPs: {str(e)}")
            raise

    def identify_significant_snps(self, p_threshold: float = 0.05) -> List[Dict]:
        """
        Identify significant trait-associated SNPs.

        Args:
            p_threshold: P-value threshold for significance

        Returns:
            List of significant SNPs
        """
        try:
            significant = self.snp_stats[self.snp_stats['chi2_p_value'] < p_threshold]
            snps = significant.to_dict('records')

            logger.info(f"Identified {len(snps)} significant trait-associated SNPs")
            return snps

        except Exception as e:
            logger.error(f"Error identifying significant SNPs: {str(e)}")
            raise

    def generate_snp_report(self) -> Dict:
        """
        Generate summary statistics for SNP analysis.

        Returns:
            Dictionary containing SNP summary statistics
        """
        try:
            stats = {
                'total_snps': len(self.snp_stats),
                'significant_snps': len(self.snp_stats[
                    self.snp_stats['chi2_p_value'] < 0.05
                ]),
                'mean_chi2': float(self.snp_stats['chi2_statistic'].mean()),
                'median_p_value': float(self.snp_stats['chi2_p_value'].median())
            }

            logger.info("Generated SNP analysis report")
            return stats

        except Exception as e:
            logger.error(f"Error generating SNP report: {str(e)}")
            raise
