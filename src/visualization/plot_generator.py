"""
Plot Generation Implementation.
"""
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

logger = logging.getLogger(__name__)

class PlotGenerator:
    """Generates visualizations for genetic analysis results."""

    def __init__(self):
        """Initialize PlotGenerator."""
        self.output_dir = Path('plots')
        self.output_dir.mkdir(exist_ok=True)

        # Set style
        plt.style.use('default')  # Use matplotlib's default style instead of seaborn
        sns.set_palette("husl")   # Still use seaborn's color palette

    def create_marker_density_heatmap(self, marker_data: pd.DataFrame,
                                    chromosomes: List[str],
                                    output_file: str) -> str:
        """
        Create heatmap of marker density across chromosomes.

        Args:
            marker_data: DataFrame with marker positions and densities
            chromosomes: List of chromosome names
            output_file: Output file path

        Returns:
            Path to saved plot
        """
        try:
            plt.figure(figsize=(12, 8))

            # Create density matrix
            density_matrix = marker_data.pivot(
                index='position_bin',
                columns='chromosome',
                values='density'
            )

            # Generate heatmap
            sns.heatmap(
                density_matrix,
                cmap='YlOrRd',
                xticklabels=chromosomes,
                cbar_kws={'label': 'Marker Density'}
            )

            plt.title('Gene Marker Density Across Chromosomes')
            plt.xlabel('Chromosome')
            plt.ylabel('Position (Mb)')

            # Save plot
            output_path = self.output_dir / output_file
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()

            logger.info(f"Created marker density heatmap: {output_path}")
            return str(output_path)

        except Exception as e:
            logger.error(f"Error creating marker density heatmap: {str(e)}")
            raise

    def plot_snp_distribution(self, snp_data: pd.DataFrame,
                            populations: List[str],
                            output_file: str) -> str:
        """
        Create SNP distribution plot across populations.

        Args:
            snp_data: DataFrame with SNP frequencies by population
            populations: List of population names
            output_file: Output file path

        Returns:
            Path to saved plot
        """
        try:
            plt.figure(figsize=(10, 6))

            # Create violin plots for SNP distribution
            sns.violinplot(
                data=snp_data,
                x='population',
                y='frequency',
                inner='box'
            )

            plt.title('SNP Frequency Distribution by Population')
            plt.xlabel('Population')
            plt.ylabel('Allele Frequency')
            plt.xticks(rotation=45)

            # Save plot
            output_path = self.output_dir / output_file
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()

            logger.info(f"Created SNP distribution plot: {output_path}")
            return str(output_path)

        except Exception as e:
            logger.error(f"Error creating SNP distribution plot: {str(e)}")
            raise

    def create_trait_distribution_map(self, trait_data: pd.DataFrame,
                                    populations: List[str],
                                    traits: List[str],
                                    output_file: str) -> str:
        """
        Create multi-population trait distribution map.

        Args:
            trait_data: DataFrame with trait values by population
            populations: List of population names
            traits: List of trait names
            output_file: Output file path

        Returns:
            Path to saved plot
        """
        try:
            n_traits = len(traits)
            fig, axes = plt.subplots(1, n_traits, figsize=(6*n_traits, 6))

            if n_traits == 1:
                axes = [axes]

            for ax, trait in zip(axes, traits):
                # Create boxplot for each trait
                sns.boxplot(
                    data=trait_data,
                    x='population',
                    y=trait,
                    ax=ax
                )

                ax.set_title(f'{trait} Distribution')
                ax.set_xlabel('Population')
                ax.set_ylabel('Value')
                ax.tick_params(axis='x', rotation=45)

            plt.tight_layout()

            # Save plot
            output_path = self.output_dir / output_file
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()


            logger.info(f"Created trait distribution map: {output_path}")
            return str(output_path)

        except Exception as e:
            logger.error(f"Error creating trait distribution map: {str(e)}")
            raise

    def plot_diversity_indices(self, diversity_data: pd.DataFrame,
                             output_file: str) -> str:
        """
        Create diversity indices comparison plot.

        Args:
            diversity_data: DataFrame with diversity indices by population
            output_file: Output file path

        Returns:
            Path to saved plot
        """
        try:
            plt.figure(figsize=(12, 6))

            # Create grouped bar plot
            indices = ['heterozygosity_observed', 'heterozygosity_expected']
            x = np.arange(len(diversity_data))
            width = 0.35

            plt.bar(x - width/2, diversity_data['heterozygosity_observed'],
                   width, label='Observed')
            plt.bar(x + width/2, diversity_data['heterozygosity_expected'],
                   width, label='Expected')

            plt.xlabel('Population')
            plt.ylabel('Heterozygosity')
            plt.title('Genetic Diversity Indices by Population')
            plt.xticks(x, diversity_data['population'], rotation=45)
            plt.legend()

            # Save plot
            output_path = self.output_dir / output_file
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()

            logger.info(f"Created diversity indices plot: {output_path}")
            return str(output_path)

        except Exception as e:
            logger.error(f"Error creating diversity indices plot: {str(e)}")
            raise
