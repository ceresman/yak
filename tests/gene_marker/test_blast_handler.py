"""
Tests for BLAST handler functionality.
"""
import os
import pytest
import pandas as pd
from pathlib import Path
from src.gene_marker.blast_handler import BlastHandler

@pytest.fixture
def blast_handler():
    """Create a BlastHandler instance for testing."""
    return BlastHandler(output_dir="tests/data/blast_output")

@pytest.fixture
def sample_alignments():
    """Create sample alignment results for testing."""
    return [
        {
            'query_id': 'gene1',
            'subject_id': 'marker1',
            'identity': 95.0,
            'alignment_length': 100,
            'e_value': 1e-50,
            'bit_score': 190.0,
            'query_start': 1,
            'query_end': 100,
            'subject_start': 1,
            'subject_end': 100
        },
        {
            'query_id': 'gene2',
            'subject_id': 'marker2',
            'identity': 85.0,
            'alignment_length': 200,
            'e_value': 1e-30,
            'bit_score': 150.0,
            'query_start': 1,
            'query_end': 200,
            'subject_start': 1,
            'subject_end': 200
        }
    ]

def test_blast_handler_init(blast_handler):
    """Test BlastHandler initialization."""
    assert isinstance(blast_handler, BlastHandler)
    assert blast_handler.output_dir.exists()

def test_calculate_similarity_scores(blast_handler, sample_alignments):
    """Test similarity score calculation."""
    df = blast_handler.calculate_similarity_scores(sample_alignments)

    assert isinstance(df, pd.DataFrame)
    assert 'similarity_score' in df.columns
    assert 'normalized_score' in df.columns
    assert len(df) == len(sample_alignments)
    assert all(0 <= score <= 1 for score in df['normalized_score'])

def test_annotate_genes(blast_handler, sample_alignments):
    """Test gene annotation."""
    similarity_df = blast_handler.calculate_similarity_scores(sample_alignments)
    annotated_df, stats = blast_handler.annotate_genes(similarity_df, threshold=0.5)

    assert isinstance(annotated_df, pd.DataFrame)
    assert isinstance(stats, dict)
    assert 'is_significant' in annotated_df.columns
    assert all(key in stats for key in ['total_alignments', 'significant_alignments',
                                      'mean_similarity', 'std_similarity'])

def test_align_sequences(blast_handler):
    """Test BLAST alignment with sample data."""
    query_file = "data/raw/markers/trait_markers.fasta"
    subject_file = "data/raw/markers/trait_markers.fasta"

    if not os.path.exists(query_file) or not os.path.exists(subject_file):
        pytest.skip("Sample FASTA files not available")

    alignments = list(blast_handler.align_sequences(query_file, subject_file))
    assert len(alignments) > 0
    assert all(isinstance(alignment, dict) for alignment in alignments)
    assert all(key in alignments[0] for key in ['query_id', 'subject_id', 'identity',
                                               'alignment_length', 'e_value'])
