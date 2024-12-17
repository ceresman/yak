"""
Tests for BLAST handler functionality.
"""
import pytest
import pandas as pd
from unittest.mock import patch, MagicMock
from src.gene_marker.blast_handler import BlastHandler

@pytest.fixture
def blast_handler(tmp_path):
    """Create a BlastHandler instance for testing."""
    output_dir = tmp_path / "blast_output"
    output_dir.mkdir(exist_ok=True)
    return BlastHandler(output_dir=str(output_dir))

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

@patch('Bio.Blast.Applications.NcbiblastnCommandline')
@patch('Bio.Blast.NCBIXML.parse')
def test_align_sequences(mock_parse, mock_blast, blast_handler, tmp_path):
    """Test BLAST alignment with mocked BLAST command."""
    # Mock BLAST command execution
    mock_cmd = MagicMock()
    mock_cmd.return_value = ('', '')
    mock_blast.return_value = mock_cmd

    # Create test files
    query_file = tmp_path / "query.fasta"
    subject_file = tmp_path / "subject.fasta"
    query_file.write_text(">seq1\nATCGATCGAT\n")
    subject_file.write_text(">ref1\nATCGATCGAT\n")

    # Set up mock BLAST record
    mock_record = MagicMock()
    mock_record.query = "seq1"

    # Set up mock alignment
    mock_alignment = MagicMock()
    mock_alignment.hit_def = "ref1"

    # Set up mock HSP
    mock_hsp = MagicMock()
    mock_hsp.identities = 95
    mock_hsp.align_length = 100
    mock_hsp.expect = 1e-50
    mock_hsp.bits = 190.0
    mock_hsp.query_start = 1
    mock_hsp.query_end = 100
    mock_hsp.sbjct_start = 1
    mock_hsp.sbjct_end = 100

    # Set up relationships
    mock_alignment.hsps = [mock_hsp]
    mock_record.alignments = [mock_alignment]

    # Configure mock parse to return our mock record
    mock_parse.return_value = iter([mock_record])

    # Run alignment
    alignments = list(blast_handler.align_sequences(str(query_file), str(subject_file)))

    # Verify results
    assert len(alignments) == 1
    assert alignments[0]['query_id'] == "seq1"
    assert alignments[0]['subject_id'] == "ref1"
    assert alignments[0]['identity'] == 95.0
