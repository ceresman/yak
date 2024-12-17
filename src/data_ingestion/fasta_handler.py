"""
FASTA file handler for processing genomic sequence data.
"""
from pathlib import Path
from typing import Generator, Optional
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class FastaHandler:
    """Handles loading and validation of FASTA files containing genomic sequences."""

    def __init__(self, cache_dir: Optional[Path] = None):
        """Initialize the FASTA handler.

        Args:
            cache_dir: Optional directory for caching processed sequences
        """
        self.cache_dir = cache_dir
        if cache_dir:
            cache_dir.mkdir(parents=True, exist_ok=True)

    def load_fasta(self, filepath: Path) -> Generator[SeqRecord, None, None]:
        """Load and validate FASTA file.

        Args:
            filepath: Path to the FASTA file

        Returns:
            Generator yielding validated SeqRecord objects

        Raises:
            FileNotFoundError: If the file doesn't exist
            ValueError: If the file is empty or invalid
        """
        try:
            if not filepath.exists():
                raise FileNotFoundError(f"FASTA file not found: {filepath}")

            logger.info(f"Loading FASTA file: {filepath}")
            sequences = SeqIO.parse(filepath, "fasta")

            # Validate file is not empty
            try:
                first_record = next(sequences)
                yield first_record
            except StopIteration:
                raise ValueError(f"Empty FASTA file: {filepath}")

            # Process remaining sequences
            for record in sequences:
                # Basic validation
                if not record.seq or len(record.seq) == 0:
                    logger.warning(f"Skipping empty sequence: {record.id}")
                    continue

                if not record.id:
                    logger.warning("Skipping sequence with no ID")
                    continue

                yield record

        except Exception as e:
            logger.error(f"Error processing FASTA file {filepath}: {str(e)}")
            raise

    def filter_contaminants(self, sequences: Generator[SeqRecord, None, None],
                          min_length: int = 100) -> Generator[SeqRecord, None, None]:
        """Filter out potential contaminant sequences.

        Args:
            sequences: Generator of SeqRecord objects
            min_length: Minimum sequence length to keep

        Returns:
            Generator yielding filtered sequences
        """
        for record in sequences:
            # Skip sequences that are too short
            if len(record.seq) < min_length:
                logger.info(f"Filtering short sequence: {record.id} (length={len(record.seq)})")
                continue

            # Add more contamination checks here as needed
            yield record
