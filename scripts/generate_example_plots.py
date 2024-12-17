"""Generate example plots using the PlotGenerator class."""
import logging
import numpy as np
import pandas as pd
from pathlib import Path
import sys
from typing import Tuple, List

# Add project root to Python path
sys.path.append(str(Path(__file__).parent.parent))

from src.visualization.plot_generator import PlotGenerator

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def generate_marker_density_data() -> Tuple[pd.DataFrame, List[str]]:
    """Generate sample marker density data."""
    chromosomes = [f'Chr{i}' for i in range(1, 21)]
    positions = np.arange(0, 100, 5)  # 20 position bins
    data = []

    for chrom in chromosomes:
        for pos in positions:
            density = np.random.poisson(lam=5)  # Random density values
            data.append({
                'chromosome': chrom,
                'position_bin': pos,
                'density': density
            })

    return pd.DataFrame(data), chromosomes

def generate_snp_data() -> Tuple[pd.DataFrame, List[str]]:
    """Generate sample SNP frequency data."""
    populations = ['Highland', 'Lowland', 'Valley', 'Mountain']
    data = []

    for pop in populations:
        # Generate 100 SNPs per population with varying frequencies
        frequencies = np.random.beta(2, 5, size=100)
        data.extend([{'population': pop, 'frequency': freq} for freq in frequencies])

    return pd.DataFrame(data), populations

def generate_trait_data() -> Tuple[pd.DataFrame, List[str], List[str]]:
    """Generate sample trait data."""
    populations = ['Highland', 'Lowland', 'Valley', 'Mountain']
    traits = ['Milk Production', 'Body Weight', 'Disease Resistance']
    data = []

    for pop in populations:
        for _ in range(30):  # 30 samples per population
            row = {'population': pop}
            row['Milk Production'] = np.random.normal(25, 5)
            row['Body Weight'] = np.random.normal(500, 50)
            row['Disease Resistance'] = np.random.beta(5, 2)
            data.append(row)

    return pd.DataFrame(data), populations, traits

def generate_diversity_data() -> pd.DataFrame:
    """Generate sample diversity data."""
    populations = ['Highland', 'Lowland', 'Valley', 'Mountain']
    data = []

    for pop in populations:
        obs_het = np.random.uniform(0.3, 0.7)
        exp_het = obs_het + np.random.uniform(-0.1, 0.1)
        data.append({
            'population': pop,
            'heterozygosity_observed': obs_het,
            'heterozygosity_expected': exp_het
        })

    return pd.DataFrame(data)

def main():
    """Generate all example plots."""
    # Initialize PlotGenerator
    plot_gen = PlotGenerator()

    try:
        # Generate marker density plot
        marker_data, chromosomes = generate_marker_density_data()
        plot_gen.create_marker_density_heatmap(
            marker_data,
            chromosomes,
            'marker_density.png'
        )

        # Generate SNP distribution plot
        snp_data, populations = generate_snp_data()
        plot_gen.plot_snp_distribution(
            snp_data,
            populations,
            'snp_distribution.png'
        )

        # Generate trait distribution plot
        trait_data, populations, traits = generate_trait_data()
        plot_gen.create_trait_distribution_map(
            trait_data,
            populations,
            traits,
            'trait_distribution.png'
        )

        # Generate diversity indices plot
        diversity_data = generate_diversity_data()
        plot_gen.plot_diversity_indices(
            diversity_data,
            'diversity_indices.png'
        )

        logger.info("All example plots generated successfully!")

    except Exception as e:
        logger.error(f"Error generating plots: {str(e)}")
        raise

if __name__ == "__main__":
    main()
