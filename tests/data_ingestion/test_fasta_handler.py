"""
Tests for FASTA file handling functionality.
"""
import os
from pathlib import Path
import pytest
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from src.data_ingestion.fasta_handler import FastaHandler

@pytest.fixture
def sample_fasta(tmp_path):
    """Create a sample FASTA file for testing."""
    fasta_file = tmp_path / "test.fasta"
    sequences = [
        SeqRecord(
            Seq("ATCGATCGATCG"),
            id="seq1",
            description="test sequence 1"
        ),
        SeqRecord(
            Seq("GCTAGCTAGCTA"),
            id="seq2",
            description="test sequence 2"
        )
    ]
    SeqIO.write(sequences, fasta_file, "fasta")
    return str(fasta_file)

@pytest.fixture
def fasta_handler():
    """Create FastaHandler instance."""
    return FastaHandler()

def test_fasta_handler_init(fasta_handler):
    """Test FastaHandler initialization."""
    assert isinstance(fasta_handler, FastaHandler)

def test_load_fasta(fasta_handler, sample_fasta):
    """Test loading FASTA file."""
    sequences = fasta_handler.load_fasta(sample_fasta)
    assert len(sequences) == 2
    assert sequences[0].id == "seq1"
    assert str(sequences[0].seq) == "ATCGATCGATCG"

def test_filter_sequences(fasta_handler, sample_fasta):
    """Test sequence filtering."""
    sequences = fasta_handler.load_fasta(sample_fasta)
    filtered = fasta_handler.filter_sequences(sequences, min_length=10)
    assert len(filtered) == 2
    filtered = fasta_handler.filter_sequences(sequences, min_length=15)
    assert len(filtered) == 0

def test_remove_ambiguous(fasta_handler):
    """Test removing ambiguous nucleotides."""
    seq = SeqRecord(Seq("ATCGNNATCG"), id="test")
    cleaned = fasta_handler.remove_ambiguous([seq])[0]
    assert str(cleaned.seq) == "ATCGATCG"

def test_error_handling(fasta_handler, tmp_path):
    """Test error handling for invalid files."""
    invalid_file = tmp_path / "nonexistent.fasta"
    with pytest.raises(FileNotFoundError):
        fasta_handler.load_fasta(str(invalid_file))
