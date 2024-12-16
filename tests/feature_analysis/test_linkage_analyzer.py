"""
Tests for linkage disequilibrium analysis functionality.
"""
import numpy as np
import pandas as pd
import pytest
from src.feature_analysis.linkage_analyzer import LinkageAnalyzer

@pytest.fixture
def sample_genotypes():
    """Create sample genotype data for testing."""
    np.random.seed(42)
    n_samples = 100
    n_markers = 5

    # Generate random genotypes (0, 1, 2 for diploid organisms)
    genotypes = np.random.choice([0, 1, 2], size=(n_samples, n_markers))
    return pd.DataFrame(genotypes, columns=[f'marker_{i}' for i in range(n_markers)])

def test_linkage_analyzer_init():
    """Test LinkageAnalyzer initialization."""
    analyzer = LinkageAnalyzer()
    assert isinstance(analyzer, LinkageAnalyzer)
    assert isinstance(analyzer.ld_results, pd.DataFrame)
    assert len(analyzer.ld_results) == 0

def test_calculate_ld(sample_genotypes):
    """Test LD calculation."""
    analyzer = LinkageAnalyzer()
    results = analyzer.calculate_ld(sample_genotypes)

    assert isinstance(results, pd.DataFrame)
    assert len(results) == (5 * 4) // 2  # n*(n-1)/2 pairs for 5 markers
    assert all(col in results.columns for col in ['marker1', 'marker2', 'D', 'D_prime', 'r2'])
    assert all(0 <= d <= 1 for d in results['D_prime'])
    assert all(0 <= r2 <= 1 for r2 in results['r2'])

def test_identify_trait_associations(sample_genotypes):
    """Test trait association identification."""
    analyzer = LinkageAnalyzer()
    analyzer.calculate_ld(sample_genotypes)

    associations = analyzer.identify_trait_associations(ld_threshold=0.5)
    assert isinstance(associations, list)
    assert all(isinstance(assoc, dict) for assoc in associations)
    assert all(assoc['D_prime'] >= 0.5 for assoc in associations)

def test_generate_ld_report(sample_genotypes):
    """Test LD report generation."""
    analyzer = LinkageAnalyzer()
    analyzer.calculate_ld(sample_genotypes)

    report = analyzer.generate_ld_report()
    assert isinstance(report, dict)
    assert all(key in report for key in ['total_pairs', 'mean_d_prime', 'mean_r2', 'high_ld_pairs'])
    assert report['total_pairs'] > 0
    assert 0 <= report['mean_d_prime'] <= 1
    assert 0 <= report['mean_r2'] <= 1
