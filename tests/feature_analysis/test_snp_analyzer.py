"""
Tests for SNP analysis functionality.
"""
import numpy as np
import pandas as pd
import pytest
from src.feature_analysis.snp_analyzer import SNPAnalyzer

@pytest.fixture
def sample_data():
    """Create sample SNP and phenotype data for testing."""
    np.random.seed(42)
    n_samples = 100
    n_snps = 5

    # Generate random SNP genotypes
    genotypes = pd.DataFrame(
        np.random.choice([0, 1, 2], size=(n_samples, n_snps)),
        columns=[f'snp_{i}' for i in range(n_snps)]
    )

    # Generate random phenotypes
    phenotypes = pd.DataFrame({
        'trait_value': np.random.choice(['case', 'control'], size=n_samples),
        'population': np.random.choice(['pop1', 'pop2'], size=n_samples)
    })

    return genotypes, phenotypes

def test_snp_analyzer_init():
    """Test SNPAnalyzer initialization."""
    analyzer = SNPAnalyzer()
    assert isinstance(analyzer, SNPAnalyzer)
    assert isinstance(analyzer.snp_stats, pd.DataFrame)
    assert len(analyzer.snp_stats) == 0

def test_evaluate_trait_snps(sample_data):
    """Test SNP evaluation."""
    genotypes, phenotypes = sample_data
    analyzer = SNPAnalyzer()

    stats = analyzer.evaluate_trait_snps(genotypes, phenotypes)
    assert isinstance(stats, pd.DataFrame)
    assert len(stats) == len(genotypes.columns)
    assert all(col in stats.columns for col in [
        'snp', 'chi2_statistic', 'chi2_p_value', 'fisher_p_value',
        'odds_ratio', 'allele_freqs'
    ])
    assert all(isinstance(freq, dict) for freq in stats['allele_freqs'])

def test_identify_significant_snps(sample_data):
    """Test significant SNP identification."""
    genotypes, phenotypes = sample_data
    analyzer = SNPAnalyzer()
    analyzer.evaluate_trait_snps(genotypes, phenotypes)

    snps = analyzer.identify_significant_snps(p_threshold=0.1)
    assert isinstance(snps, list)
    assert all(isinstance(snp, dict) for snp in snps)
    assert all(snp['chi2_p_value'] < 0.1 for snp in snps)

def test_generate_snp_report(sample_data):
    """Test SNP report generation."""
    genotypes, phenotypes = sample_data
    analyzer = SNPAnalyzer()
    analyzer.evaluate_trait_snps(genotypes, phenotypes)

    report = analyzer.generate_snp_report()
    assert isinstance(report, dict)
    assert all(key in report for key in [
        'total_snps', 'significant_snps', 'mean_chi2', 'median_p_value'
    ])
    assert report['total_snps'] > 0
    assert report['significant_snps'] >= 0
    assert report['mean_chi2'] > 0
