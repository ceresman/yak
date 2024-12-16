"""
FASTA file handling functionality.
"""
from pathlib import Path
from typing import List
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class FastaHandler:
    """Handles FASTA file operations for genomic data."""

    def __init__(self):
        """Initialize FastaHandler."""
        self.logger = logger

    def load_fasta(self, file_path: str) -> List[SeqRecord]:
        """
        Load sequences from a FASTA file.

        Args:
            file_path: Path to the FASTA file

        Returns:
            List of SeqRecord objects
        """
        try:
            self.logger.info(f"Loading FASTA file: {file_path}")
            sequences = list(SeqIO.parse(file_path, "fasta"))
            self.logger.info(f"Successfully loaded {len(sequences)} sequences")
            return sequences
        except FileNotFoundError:
            self.logger.error(f"FASTA file not found: {file_path}")
            raise
        except Exception as e:
            self.logger.error(f"Error loading FASTA file: {str(e)}")
            raise

    def filter_sequences(self, sequences: List[SeqRecord], min_length: int = 0) -> List[SeqRecord]:
        """
        Filter sequences based on length.

        Args:
            sequences: List of SeqRecord objects
            min_length: Minimum sequence length to keep

        Returns:
            Filtered list of SeqRecord objects
        """
        try:
            filtered = [seq for seq in sequences if len(seq.seq) >= min_length]
            self.logger.info(f"Filtered {len(sequences) - len(filtered)} sequences below length {min_length}")
            return filtered
        except Exception as e:
            self.logger.error(f"Error filtering sequences: {str(e)}")
            raise

    def remove_ambiguous(self, sequences: List[SeqRecord]) -> List[SeqRecord]:
        """
        Remove ambiguous nucleotides from sequences.

        Args:
            sequences: List of SeqRecord objects

        Returns:
            List of SeqRecord objects with ambiguous nucleotides removed
        """
        try:
            cleaned = []
            for seq in sequences:
                # Remove N's and other ambiguous nucleotides
                clean_seq = str(seq.seq).replace('N', '').replace('n', '')
                cleaned.append(SeqRecord(Seq(clean_seq), id=seq.id, description=seq.description))

            self.logger.info(f"Cleaned ambiguous nucleotides from {len(sequences)} sequences")
            return cleaned
        except Exception as e:
            self.logger.error(f"Error removing ambiguous nucleotides: {str(e)}")
            raise
