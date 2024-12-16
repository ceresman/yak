"""
Tests for plot generation functionality.
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from src.visualization.plot_generator import PlotGenerator

@pytest.fixture
def sample_data():
    """Create sample data for testing plots."""
    np.random.seed(42)

    # Marker density data
    marker_data = pd.DataFrame({
        'chromosome': np.repeat(['chr1', 'chr2', 'chr3'], 10),
        'position_bin': np.tile(range(10), 3),
        'density': np.random.rand(30)
    })

    # SNP frequency data
    snp_data = pd.DataFrame({
        'population': np.repeat(['pop1', 'pop2', 'pop3'], 100),
        'frequency': np.random.rand(300)
    })

    # Trait distribution data
    trait_data = pd.DataFrame({
        'population': np.repeat(['pop1', 'pop2', 'pop3'], 50),
        'trait1': np.random.normal(0, 1, 150),
        'trait2': np.random.normal(0, 1, 150)
    })

    # Diversity data
    diversity_data = pd.DataFrame({
        'population': ['pop1', 'pop2', 'pop3'],
        'heterozygosity_observed': np.random.rand(3),
        'heterozygosity_expected': np.random.rand(3)
    })

    return marker_data, snp_data, trait_data, diversity_data

def test_plot_generator_init():
    """Test PlotGenerator initialization."""
    generator = PlotGenerator()
    assert isinstance(generator, PlotGenerator)
    assert generator.output_dir.exists()

def test_create_marker_density_heatmap(sample_data):
    """Test marker density heatmap creation."""
    marker_data, _, _, _ = sample_data
    generator = PlotGenerator()

    output_file = 'test_heatmap.png'
    plot_path = generator.create_marker_density_heatmap(
        marker_data,
        ['chr1', 'chr2', 'chr3'],
        output_file
    )

    assert os.path.exists(plot_path)
    assert plot_path.endswith('.png')

def test_plot_snp_distribution(sample_data):
    """Test SNP distribution plot creation."""
    _, snp_data, _, _ = sample_data
    generator = PlotGenerator()

    output_file = 'test_snp_dist.png'
    plot_path = generator.plot_snp_distribution(
        snp_data,
        ['pop1', 'pop2', 'pop3'],
        output_file
    )

    assert os.path.exists(plot_path)
    assert plot_path.endswith('.png')

def test_create_trait_distribution_map(sample_data):
    """Test trait distribution map creation."""
    _, _, trait_data, _ = sample_data
    generator = PlotGenerator()

    output_file = 'test_trait_dist.png'
    plot_path = generator.create_trait_distribution_map(
        trait_data,
        ['pop1', 'pop2', 'pop3'],
        ['trait1', 'trait2'],
        output_file
    )

    assert os.path.exists(plot_path)
    assert plot_path.endswith('.png')

def test_plot_diversity_indices(sample_data):
    """Test diversity indices plot creation."""
    _, _, _, diversity_data = sample_data
    generator = PlotGenerator()

    output_file = 'test_diversity.png'
    plot_path = generator.plot_diversity_indices(
        diversity_data,
        output_file
    )

    assert os.path.exists(plot_path)
    assert plot_path.endswith('.png')

@pytest.fixture(autouse=True)
def cleanup():
    """Clean up test plots after each test."""
    yield
    plot_dir = Path('plots')
    if plot_dir.exists():
        for file in plot_dir.glob('*.png'):
            file.unlink()
        plot_dir.rmdir()
