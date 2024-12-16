"""
Tests for BLAST handler functionality.
"""
import os
import pytest
import pandas as pd
from unittest.mock import patch, MagicMock
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
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

@pytest.fixture
def mock_blast_xml(tmp_path):
    """Create mock BLAST XML output."""
    xml_content = """<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_hits>
        <Hit>
          <Hit_def>marker1</Hit_def>
          <Hit_hsps>
            <Hsp>
              <Hsp_bit-score>190.0</Hsp-bit-score>
              <Hsp_score>95</Hsp_score>
              <Hsp_expect>1e-50</Hsp_expect>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>100</Hsp_query-to>
              <Hsp_hit-from>1</Hsp_hit-from>
              <Hsp_hit-to>100</Hsp_hit-to>
              <Hsp_identity>95</Hsp_identity>
              <Hsp_align-len>100</Hsp_align-len>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>
"""
    xml_file = tmp_path / "blast_results.xml"
    xml_file.write_text(xml_content)
    return str(xml_file)

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
def test_align_sequences(mock_blast, blast_handler, mock_blast_xml, tmp_path):
    """Test BLAST alignment with mocked BLAST command."""
    # Mock BLAST command execution
    mock_cmd = MagicMock()
    mock_cmd.return_value = ('', '')
    mock_blast.return_value = mock_cmd

    # Create test files
    query_file = tmp_path / "query.fasta"
    subject_file = tmp_path / "subject.fasta"
    query_file.touch()
    subject_file.touch()

    # Run alignment with mocked BLAST
    alignments = list(blast_handler.align_sequences(str(query_file), str(subject_file)))

    assert len(alignments) > 0
    assert all(isinstance(alignment, dict) for alignment in alignments)
    assert all(key in alignments[0] for key in ['query_id', 'subject_id', 'identity',
                                              'alignment_length', 'e_value'])
