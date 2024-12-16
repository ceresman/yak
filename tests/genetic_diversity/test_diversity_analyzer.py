"""
Tests for genetic diversity analysis functionality.
"""
import numpy as np
import pandas as pd
import pytest
from src.genetic_diversity.diversity_analyzer import DiversityAnalyzer

@pytest.fixture
def sample_data():
    """Create sample genotype and population data for testing."""
    np.random.seed(42)
    n_samples = 100
    n_markers = 5

    # Generate random genotypes
    genotypes = pd.DataFrame(
        np.random.choice([0, 1, 2], size=(n_samples, n_markers)),
        columns=[f'marker_{i}' for i in range(n_markers)]
    )

    # Generate population labels
    populations = pd.Series(
        np.random.choice(['pop1', 'pop2', 'pop3'], size=n_samples)
    )

    return genotypes, populations

def test_diversity_analyzer_init():
    """Test DiversityAnalyzer initialization."""
    analyzer = DiversityAnalyzer()
    assert isinstance(analyzer, DiversityAnalyzer)
    assert isinstance(analyzer.diversity_stats, pd.DataFrame)
    assert len(analyzer.diversity_stats) == 0

def test_calculate_diversity_indices(sample_data):
    """Test diversity indices calculation."""
    genotypes, populations = sample_data
    analyzer = DiversityAnalyzer()

    stats = analyzer.calculate_diversity_indices(genotypes, populations)
    assert isinstance(stats, pd.DataFrame)
    assert len(stats) == len(populations.unique())
    assert all(col in stats.columns for col in [
        'population', 'sample_size', 'mean_heterozygosity_expected',
        'mean_heterozygosity_observed', 'effective_size', 'polymorphic_loci'
    ])
    assert all(0 <= stats['mean_heterozygosity_expected'])
    assert all(stats['mean_heterozygosity_expected'] <= 1)

def test_analyze_population_structure(sample_data):
    """Test population structure analysis."""
    genotypes, populations = sample_data
    analyzer = DiversityAnalyzer()

    structure = analyzer.analyze_population_structure(genotypes, populations)
    assert isinstance(structure, dict)
    assert all(key in structure for key in [
        'global_fst', 'pairwise_fst', 'population_count', 'mean_pairwise_fst'
    ])
    assert 0 <= structure['global_fst'] <= 1
    assert all(0 <= fst <= 1 for fst in structure['pairwise_fst'].values())

def test_detect_bottlenecks(sample_data):
    """Test bottleneck detection."""
    genotypes, populations = sample_data
    analyzer = DiversityAnalyzer()

    bottlenecks = analyzer.detect_bottlenecks(genotypes, populations)
    assert isinstance(bottlenecks, dict)
    assert len(bottlenecks) == len(populations.unique())
    for pop_stats in bottlenecks.values():
        assert all(key in pop_stats for key in [
            'heterozygosity_excess', 'wilcoxon_p_value', 'low_freq_alleles'
        ])
        assert isinstance(pop_stats['wilcoxon_p_value'], float)
        assert 0 <= pop_stats['low_freq_alleles'] <= 1
