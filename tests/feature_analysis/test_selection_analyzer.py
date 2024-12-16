"""
Tests for selection marker analysis functionality.
"""
import numpy as np
import pandas as pd
import pytest
from src.feature_analysis.selection_analyzer import SelectionAnalyzer

@pytest.fixture
def sample_data():
    """Create sample genotype and phenotype data for testing."""
    np.random.seed(42)
    n_samples = 100
    n_markers = 5

    # Generate random genotypes
    genotypes = pd.DataFrame(
        np.random.choice([0, 1, 2], size=(n_samples, n_markers)),
        columns=[f'marker_{i}' for i in range(n_markers)]
    )

    # Generate random phenotypes
    phenotypes = pd.DataFrame({
        'trait': np.random.choice(['A', 'B'], size=n_samples),
        'population': np.random.choice(['pop1', 'pop2'], size=n_samples)
    })

    return genotypes, phenotypes

def test_selection_analyzer_init():
    """Test SelectionAnalyzer initialization."""
    analyzer = SelectionAnalyzer()
    assert isinstance(analyzer, SelectionAnalyzer)
    assert isinstance(analyzer.selection_stats, pd.DataFrame)
    assert len(analyzer.selection_stats) == 0

def test_calculate_selection_statistics(sample_data):
    """Test selection statistics calculation."""
    genotypes, phenotypes = sample_data
    analyzer = SelectionAnalyzer()

    stats = analyzer.calculate_selection_statistics(genotypes, phenotypes)
    assert isinstance(stats, pd.DataFrame)
    assert len(stats) == len(genotypes.columns)
    assert all(col in stats.columns for col in [
        'marker', 'allele_frequency', 'heterozygosity', 'fst',
        'chi2_statistic', 'p_value'
    ])
    assert all(0 <= freq <= 1 for freq in stats['allele_frequency'])
    assert all(0 <= het <= 1 for het in stats['heterozygosity'])
    assert all(0 <= fst <= 1 for fst in stats['fst'])

def test_identify_selection_regions(sample_data):
    """Test selection region identification."""
    genotypes, phenotypes = sample_data
    analyzer = SelectionAnalyzer()
    analyzer.calculate_selection_statistics(genotypes, phenotypes)

    regions = analyzer.identify_selection_regions(p_threshold=0.1)
    assert isinstance(regions, list)
    assert all(isinstance(region, dict) for region in regions)
    assert all(region['p_value'] < 0.1 for region in regions)

def test_generate_selection_report(sample_data):
    """Test selection report generation."""
    genotypes, phenotypes = sample_data
    analyzer = SelectionAnalyzer()
    analyzer.calculate_selection_statistics(genotypes, phenotypes)

    report = analyzer.generate_selection_report()
    assert isinstance(report, dict)
    assert all(key in report for key in [
        'total_markers', 'mean_heterozygosity', 'mean_fst',
        'significant_markers'
    ])
    assert report['total_markers'] > 0
    assert 0 <= report['mean_heterozygosity'] <= 1
    assert 0 <= report['mean_fst'] <= 1
