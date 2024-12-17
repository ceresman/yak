"""
Selection Marker Analysis Module for identifying regions under selective pressure.
"""
import logging
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from scipy import stats

logger = logging.getLogger(__name__)

class SelectionAnalyzer:
    """Analyzes genetic markers for selection and polymorphism patterns."""

    def __init__(self):
        """Initialize SelectionAnalyzer."""
        self.selection_stats = pd.DataFrame()

    def calculate_selection_statistics(self, genotypes: pd.DataFrame,
                                    phenotypes: pd.DataFrame) -> pd.DataFrame:
        """
        Calculate selection statistics for genetic markers.

        Args:
            genotypes: DataFrame with marker genotypes
            phenotypes: DataFrame with trait phenotypes

        Returns:
            DataFrame with selection statistics
        """
        try:
            selection_data = []

            for marker in genotypes.columns:
                # Calculate allele frequencies
                allele_freq = np.clip(genotypes[marker].mean() / 2, 0, 1)

                # Calculate heterozygosity
                het = np.clip(2 * allele_freq * (1 - allele_freq), 0, 1)

                # Calculate FST-like statistic for each population
                populations = phenotypes['population'].unique()
                total_var = np.var(genotypes[marker])

                # Handle zero variance case
                if total_var == 0:
                    fst = 0
                else:
                    within_vars = []
                    for pop in populations:
                        pop_mask = phenotypes['population'] == pop
                        pop_var = np.var(genotypes.loc[pop_mask, marker])
                        within_vars.append(pop_var)

                    within_var = np.mean(within_vars)
                    fst = np.clip((total_var - within_var) / total_var, 0, 1)

                # Chi-square test for selection
                # Create 2x2 contingency table by binning trait values
                binary_genotypes = (genotypes[marker] > 0).astype(int)
                binary_traits = (phenotypes['trait_value'] > phenotypes['trait_value'].median()).astype(int)
                contingency = pd.crosstab(binary_genotypes, binary_traits)

                chi2, pval = stats.chi2_contingency(contingency)[:2]

                selection_data.append({
                    'marker': marker,
                    'allele_frequency': allele_freq,
                    'heterozygosity': het,
                    'fst': fst,
                    'chi2_statistic': chi2,
                    'p_value': pval
                })

            self.selection_stats = pd.DataFrame(selection_data)
            logger.info(f"Calculated selection statistics for {len(selection_data)} markers")
            return self.selection_stats

        except Exception as e:
            logger.error(f"Error calculating selection statistics: {str(e)}")
            raise

    def identify_selection_regions(self, p_threshold: float = 0.05) -> List[Dict]:
        """
        Identify regions under selective pressure.

        Args:
            p_threshold: P-value threshold for significance

        Returns:
            List of markers under selection
        """
        try:
            selected = self.selection_stats[self.selection_stats['p_value'] < p_threshold]
            regions = selected.to_dict('records')

            logger.info(f"Identified {len(regions)} regions under selection")
            return regions

        except Exception as e:
            logger.error(f"Error identifying selection regions: {str(e)}")
            raise

    def generate_selection_report(self) -> Dict:
        """
        Generate summary statistics for selection analysis.

        Returns:
            Dictionary containing selection summary statistics
        """
        try:
            stats = {
                'total_markers': len(self.selection_stats),
                'mean_heterozygosity': float(np.clip(self.selection_stats['heterozygosity'].mean(), 0, 1)),
                'mean_fst': float(np.clip(self.selection_stats['fst'].mean(), 0, 1)),
                'significant_markers': len(self.selection_stats[
                    self.selection_stats['p_value'] < 0.05
                ])
            }

            logger.info("Generated selection analysis report")
            return stats

        except Exception as e:
            logger.error(f"Error generating selection report: {str(e)}")
            raise
