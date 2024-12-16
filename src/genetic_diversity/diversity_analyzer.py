"""
Genetic Diversity Analysis Implementation.
"""
import logging
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from scipy import stats

logger = logging.getLogger(__name__)

class DiversityAnalyzer:
    """Analyzes genetic diversity within and between populations."""

    def __init__(self):
        """Initialize DiversityAnalyzer."""
        self.diversity_stats = pd.DataFrame()

    def calculate_diversity_indices(self, genotypes: pd.DataFrame,
                                 populations: pd.Series) -> pd.DataFrame:
        """
        Calculate diversity indices for populations.

        Args:
            genotypes: DataFrame with marker genotypes
            populations: Series with population labels

        Returns:
            DataFrame with diversity statistics
        """
        try:
            diversity_data = []

            for pop in populations.unique():
                pop_mask = populations == pop
                pop_genotypes = genotypes.loc[pop_mask]

                # Calculate allele frequencies
                allele_freqs = pop_genotypes.mean() / 2

                # Calculate Nei's diversity index (He)
                he = 2 * allele_freqs * (1 - allele_freqs)
                mean_he = float(np.mean(he))

                # Calculate observed heterozygosity (Ho)
                ho = (pop_genotypes == 1).mean()
                mean_ho = float(np.mean(ho))

                # Calculate effective number of alleles (Ne)
                ne = 1 / (1 - mean_he) if mean_he < 1 else float('inf')

                diversity_data.append({
                    'population': pop,
                    'sample_size': len(pop_genotypes),
                    'mean_heterozygosity_expected': mean_he,
                    'mean_heterozygosity_observed': mean_ho,
                    'effective_size': ne,
                    'polymorphic_loci': np.sum(allele_freqs > 0) / len(allele_freqs)
                })

            self.diversity_stats = pd.DataFrame(diversity_data)
            logger.info(f"Calculated diversity indices for {len(diversity_data)} populations")
            return self.diversity_stats

        except Exception as e:
            logger.error(f"Error calculating diversity indices: {str(e)}")
            raise

    def analyze_population_structure(self, genotypes: pd.DataFrame,
                                  populations: pd.Series) -> Dict:
        """
        Analyze population structure and differentiation.

        Args:
            genotypes: DataFrame with marker genotypes
            populations: Series with population labels

        Returns:
            Dictionary with population structure statistics
        """
        try:
            # Calculate overall FST
            total_var = np.var(genotypes.values)
            within_vars = []

            for pop in populations.unique():
                pop_mask = populations == pop
                pop_var = np.var(genotypes.loc[pop_mask].values)
                within_vars.append(pop_var)

            within_var = np.mean(within_vars)
            fst = np.clip((total_var - within_var) / total_var, 0, 1)

            # Calculate pairwise FST
            pops = populations.unique()
            pairwise_fst = {}

            for i, pop1 in enumerate(pops):
                for pop2 in pops[i+1:]:
                    mask1 = populations == pop1
                    mask2 = populations == pop2
                    var1 = np.var(genotypes.loc[mask1].values)
                    var2 = np.var(genotypes.loc[mask2].values)
                    total = np.var(genotypes.loc[mask1 | mask2].values)
                    pair_fst = np.clip((total - np.mean([var1, var2])) / total, 0, 1)
                    pairwise_fst[f"{pop1}-{pop2}"] = float(pair_fst)

            structure_stats = {
                'global_fst': float(fst),
                'pairwise_fst': pairwise_fst,
                'population_count': len(pops),
                'mean_pairwise_fst': float(np.mean(list(pairwise_fst.values())))
            }

            logger.info("Analyzed population structure")
            return structure_stats

        except Exception as e:
            logger.error(f"Error analyzing population structure: {str(e)}")
            raise

    def detect_bottlenecks(self, genotypes: pd.DataFrame,
                         populations: pd.Series) -> Dict:
        """
        Detect potential genetic bottlenecks in populations.

        Args:
            genotypes: DataFrame with marker genotypes
            populations: Series with population labels

        Returns:
            Dictionary with bottleneck statistics
        """
        try:
            bottleneck_stats = {}

            for pop in populations.unique():
                pop_mask = populations == pop
                pop_genotypes = genotypes.loc[pop_mask]

                # Calculate allele frequencies
                allele_freqs = pop_genotypes.mean() / 2

                # Calculate heterozygosity excess/deficiency
                he = 2 * allele_freqs * (1 - allele_freqs)
                ho = (pop_genotypes == 1).mean()

                # Test for heterozygosity excess (sign of recent bottleneck)
                _, pval = stats.wilcoxon(ho - he)

                bottleneck_stats[pop] = {
                    'heterozygosity_excess': float(np.mean(ho - he)),
                    'wilcoxon_p_value': float(pval),
                    'low_freq_alleles': float(np.mean(allele_freqs < 0.1))
                }

            logger.info("Completed bottleneck detection analysis")
            return bottleneck_stats

        except Exception as e:
            logger.error(f"Error detecting bottlenecks: {str(e)}")
            raise
