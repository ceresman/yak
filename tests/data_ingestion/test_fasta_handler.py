"""
Tests for FASTA file handler.
"""
import pytest
from pathlib import Path
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from src.data_ingestion.fasta_handler import FastaHandler

@pytest.fixture
def sample_fasta(tmp_path):
    """Create a sample FASTA file for testing."""
    fasta_file = tmp_path / "test.fasta"
    with open(fasta_file, 'w') as f:
        f.write(">seq1\nATCG\n>seq2\nGCTA\n>empty\n\n>short\nA")
    return fasta_file

@pytest.fixture
def fasta_handler():
    """Create a FastaHandler instance."""
    return FastaHandler()

def test_load_fasta_valid_file(fasta_handler, sample_fasta):
    """Test loading a valid FASTA file."""
    sequences = list(fasta_handler.load_fasta(sample_fasta))
    assert len(sequences) == 3  # empty sequence is filtered
    assert sequences[0].id == "seq1"
    assert str(sequences[0].seq) == "ATCG"

def test_load_fasta_missing_file(fasta_handler):
    """Test loading a non-existent file."""
    with pytest.raises(FileNotFoundError):
        list(fasta_handler.load_fasta(Path("nonexistent.fasta")))

def test_filter_contaminants(fasta_handler):
    """Test filtering contaminant sequences."""
    sequences = [
        SeqRecord(Seq("ATCG" * 25), id="long_seq"),
        SeqRecord(Seq("AT"), id="short_seq")
    ]
    filtered = list(fasta_handler.filter_contaminants(iter(sequences), min_length=100))
    assert len(filtered) == 1
    assert filtered[0].id == "long_seq"
