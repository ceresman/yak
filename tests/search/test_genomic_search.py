"""
Tests for genomic search functionality.
"""
import pytest
import pandas as pd
import numpy as np
from src.search.genomic_search import GenomicSearch

@pytest.fixture
def sample_markers():
    """Create sample marker data for testing."""
    return pd.DataFrame({
        'marker_id': ['rs1', 'rs2', 'rs3', 'rs4'],
        'chromosome': ['chr1', 'chr1', 'chr2', 'chr2'],
        'position': [100, 200, 150, 300],
        'trait': ['milk', 'milk', 'size', 'disease'],
        'effect_size': [0.5, 0.3, 0.8, 0.4]
    })

def test_genomic_search_init():
    """Test GenomicSearch initialization."""
    search = GenomicSearch()
    assert isinstance(search, GenomicSearch)
    assert isinstance(search.marker_index, dict)
    assert isinstance(search.region_index, dict)

def test_index_markers(sample_markers):
    """Test marker indexing."""
    search = GenomicSearch()
    search.index_markers(sample_markers)
    assert len(search.marker_index) == 4
    assert 'rs1' in search.marker_index
    assert search.marker_index['rs1']['trait'] == 'milk'

def test_search_by_marker(sample_markers):
    """Test marker search."""
    search = GenomicSearch()
    search.index_markers(sample_markers)

    result = search.search_by_marker('rs1')
    assert result is not None
    assert result['chromosome'] == 'chr1'
    assert result['trait'] == 'milk'

    result = search.search_by_marker('nonexistent')
    assert result is None

def test_search_by_region(sample_markers):
    """Test region search."""
    search = GenomicSearch()
    search.index_markers(sample_markers)

    results = search.search_by_region('chr1', 50, 150)
    assert len(results) == 1
    assert results[0]['marker_id'] == 'rs1'

    results = search.search_by_region('chr2', 0, 1000)
    assert len(results) == 2
    assert all(r['chromosome'] == 'chr2' for r in results)

def test_search_by_trait(sample_markers):
    """Test trait search."""
    search = GenomicSearch()
    search.index_markers(sample_markers)

    results = search.search_by_trait('milk')
    assert len(results) == 2
    assert all(r['trait'] == 'milk' for r in results)

    results = search.search_by_trait('milk', min_effect=0.4)
    assert len(results) == 1
    assert results[0]['marker_id'] == 'rs1'
